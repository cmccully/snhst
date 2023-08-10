#!/usr/bin/env python
import os
from glob import glob
import shutil
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import table
from astropy import nddata
import stwcs

from snhst import drizzle, dolphot, wcs, fits_utils, utils, parameters
import reproject
import sep

import logging
logging.basicConfig(level=logging.INFO)


def run(**options):
    """Top level pipeline function to align and stack images and measure photometry.

    Inputs are flexible, but keywords can include:
     - ra, dec: coordinates of the target in decimal degrees (default = center of ground_reference)
     - ground_reference: image filename for coordinate references
     - sep: dictionary of options to sep (see sep manual)
     - dolphot: dictionary of general dolphot parameters (see dolphot manual)
     - dolphot_img: dictionary of image-specific dolphot parameters (see dolphot manual)
     - calcsky: dictionary of parameters to calcsky (r_in, r_out, step, sigma_low, sigma_high)
     - fakelist: dictionary of parameters to fakelist (filter1, filter2, filter1_min, filter1_max,
                                                       color_min, color_max, nstar)"""

    parameters.set_default_parameters(options, parameters.global_defaults)

    # Unzip images
    if glob('*.fits.gz'):
        os.system('gunzip *.fits.gz')

    # Download the reference images
    download_reference_images()

    # Copy the raw files into the raw directory
    copy_raw_data()

    # output coords from header of ground reference, if not given
    if (options['ra'] is None or options['dec'] is None) and options['ground_reference'] is not None:
        options['ra'], options['dec'] = get_target_coordinates(options['ground_reference'])

    raw_images = get_raw_image_filenames('raw')

    # Make sure the position of interest is in the raw frame
    if options['ra'] is not None and options['dec'] is not None:
        remove_images_without_object(options['ra'], options['dec'], raw_images)

    # sort the data into instrument/visit/filter
    visit_meta_data = sort_raw_data(raw_images)
    hst_reference = get_default_reference_hst(visit_meta_data)
    reference_path = os.path.join(hst_reference['instrument'], hst_reference['visit'], hst_reference['filter'])
    images = get_raw_image_filenames(reference_path)
    hst_template_path = os.path.join(reference_path, 'hst_template')
    hst_template_images = [utils.copy_if_not_exists(image, hst_template_path) for image in images]

    # Make the overall HST template image
    if options['ground_reference_catalog'] is not None or options['ground_reference'] is not None:

        match_template_path = os.path.join(reference_path, 'match_ground')
        match_template_images = [utils.copy_if_not_exists(image, match_template_path) for image in images]
        drizzled_template_filename = drizzle.drizzle(match_template_images, 'match_ground', options)

        if options['use_sep']:
            data, mask, drizzled_template_header = get_masked_data_for_sep(drizzled_template_filename)
            drizzled_template_catalog = make_sep_catalog(data, drizzled_template_header, options, mask=mask)
        else:
            dolphot_path = os.path.join(match_template_path, 'dolphot')
            drizzled_template_catalog = dolphot.dolphot(drizzled_template_filename, images, dolphot_path, options)

        if options['ground_reference_catalog'] is None:
            reference_hdu = fits_utils.get_data_extension(options['ground_reference'])
            # Reproject ground image onto drizzled HST frame
            drizzled_template_header = fits.getheader(drizzled_template_filename, extname='SCI')
            ground_reference_data, _ = reproject.reproject_interp(reference_hdu, drizzled_template_header)
            # Run sep on the ground image to make the ground catalog
            ground_reference_catalog = make_sep_catalog(ground_reference_data, drizzled_template_header, options)
        else:
            ground_reference_catalog = table.Table.read(options['ground_reference_catalog'], format='ascii')

        # Calculate the offset between the frames
        # (both catalogs are projected to the WCS of the drizzled template)
        wcs.offset_to_match(hst_template_images, drizzled_template_catalog, ground_reference_catalog, max_offset=2.)

    hst_reference_filename = drizzle.drizzle(hst_template_images, 'hst_template', options)
    if options['use_sep']:
        data, mask, hst_reference_header = get_masked_data_for_sep(hst_reference_filename)
        hst_reference_catalog = make_sep_catalog(data, hst_reference_header, options, mask=mask)
    else:
        dolphot_path = os.path.join(hst_template_path, 'dolphot')
        hst_reference_catalog = dolphot.dolphot(hst_reference_filename, images, dolphot_path, options)

    # drizzle each visit+filter separately
    offset_images_by_filter = {filt: [] for filt in visit_meta_data['filter']}
    for visit in visit_meta_data:
        visit_path = os.path.join(visit['instrument'], visit['visit'], visit['filter'])
        images = get_raw_image_filenames(visit_path)
        match_template_path = os.path.join(visit_path, 'match_template')
        match_template_images = [utils.copy_if_not_exists(image, match_template_path) for image in images]

        visit_template_filename = drizzle.drizzle(match_template_images, 'match_template', options)
        if options['use_sep']:
            data, mask, visit_template_header = get_masked_data_for_sep(visit_template_filename)
            visit_template_catalog = make_sep_catalog(data, visit_template_header, options, mask=mask)
        else:
            dolphot_path = os.path.join(match_template_path, 'dolphot')
            visit_template_catalog = dolphot.dolphot(visit_template_filename, images, dolphot_path, options)

        hst_template_path = os.path.join(visit_path, 'visit_template')
        hst_template_images = [utils.copy_if_not_exists(image, hst_template_path) for image in images]

        wcs.offset_to_match(hst_template_images, visit_template_catalog, hst_reference_catalog, max_offset=2.)
        if options['make_visit_templates']:
            visit_template_filename = drizzle.drizzle(hst_template_images, '_'.join(visit['visit', 'filter']), options)
            utils.copy_if_not_exists(visit_template_filename, 'final')

        offset_images_by_filter[visit['filter']] += hst_template_images

    # drizzle all visits together, but keep filters separate
    for filt in set(visit_meta_data['filter']):
        images = offset_images_by_filter[filt]
        hst_template_path = os.path.join(filt, 'filter_template')
        hst_template_images = [utils.copy_if_not_exists(image, hst_template_path) for image in images]

        filter_template_filename = drizzle.drizzle(hst_template_images, filt, options)
        utils.copy_if_not_exists(filter_template_filename, 'final')

    if not options['use_sep']:
        # Run dolphot on all the images together (forced photometry on same positions)
        dolphot.dolphot(hst_reference_filename, sum(offset_images_by_filter.values(), []), 'dolphot', options)
        utils.copy_if_not_exists('dolphot/hst.cat', 'final')

        # Run dolphot again with fake stars
        dolphot.add_fake_stars(visit_meta_data, 'dolphot', options['fakelist'])
        utils.copy_if_not_exists('dolphot/hst_fake.cat', 'final')


def make_sep_catalog(data, header, options, mask=None, min_sep=10., do_bgsub=False):

    try:
        bkg = sep.Background(data, mask, bw=32, bh=32, fw=3, fh=3)
    except ValueError:
        data = data.byteswap().newbyteorder()
        bkg = sep.Background(data, mask, bw=32, bh=32, fw=3, fh=3)

    if do_bgsub:
        error = np.sqrt(data)
        data_bgsub = data - bkg
    else:
        error = bkg.globalrms
        data_bgsub = data
    sources = sep.extract(data_bgsub, err=error, mask=mask, **options['sep'])

    dists = ((sources['x'] - sources['x'][:, np.newaxis])**2 +
             (sources['y'] - sources['y'][:, np.newaxis])**2)**0.5
    closest = np.partition(dists, 1)[:, 1]
    sources = sources[closest > min_sep]

    t = table.Table(sources)
    kronrad, krflag = sep.kron_radius(data_bgsub, sources['x'], sources['y'],
                                      sources['a'], sources['b'],
                                      sources['theta'], 6.0)

    flux, fluxerr, flag = sep.sum_ellipse(data_bgsub, sources['x'], sources['y'],
                                          sources['a'], sources['b'],
                                          np.pi / 2.0, 2.5 * kronrad,
                                          subpix=1, err=error)

    t['mag'] = -2.5 * np.log10(flux)
    t['magerr'] = np.log(10) / 2.5 * fluxerr / flux
    t['ra'], t['dec'] = WCS(header).all_pix2world(t['x'], t['y'], 0)

    t = t['x', 'y', 'mag', 'magerr', 'ra', 'dec']

    return t


def get_masked_data_for_sep(filename, bin_factor=32, thresh=0.7):
    hdulist = fits.open(filename)
    data = hdulist['SCI'].data
    weight = hdulist['WHT'].data
    header = hdulist['SCI'].header

    binned = nddata.block_reduce(weight, bin_factor)
    debinned = nddata.block_replicate(binned, bin_factor)

    i0 = (weight.shape[0] - debinned.shape[0]) // 2
    i1 = (debinned.shape[0] - weight.shape[0]) // 2
    j0 = (weight.shape[1] - debinned.shape[1]) // 2
    j1 = (debinned.shape[1] - weight.shape[1]) // 2

    padded = np.pad(debinned, ((i0, -i1), (j0, -j1)), 'constant')
    mask = padded < thresh * np.max(weight)

    return data, mask, header


def remove_images_without_object(ra, dec, images):
    utils.mkdir('unused')
    for image in images:
        if not image_includes_coordinate(image, ra, dec):
            shutil.move(image, 'unused/')


def image_includes_coordinate(image, ra, dec):
    image_hdu = fits.open(image)
    coordinate_in_image = False
    for i, hdu in enumerate(image_hdu):
        if hdu.header.get('EXTNAME') == 'SCI':
            try:
                image_wcs = WCS(hdu)
            except:
                image_wcs = stwcs.wcsutil.HSTWCS(image_hdu, i)
            # Use zero indexed here for the pixel coordinates
            coordinate_in_pixels = image_wcs.all_world2pix([[ra, dec]], 0)
            coordinate_not_negative = np.all(coordinate_in_pixels > 0)
            coordinate_not_outside = np.all(coordinate_in_pixels < np.array([hdu.header['NAXIS2'],
                                                                             hdu.header['NAXIS1']]))
            if coordinate_not_negative and coordinate_not_outside:
                coordinate_in_image = True
                break

    image_hdu.close()
    return coordinate_in_image


def copy_raw_data():
    fs = glob('*.fits')
    for f in fs:
        utils.copy_if_not_exists(f, 'raw/')


def download_reference_images():
    os.environ['CRDS_SERVER_URL'] = 'https://hst-crds.stsci.edu'
    os.environ['CRDS_PATH'] = os.getcwd() + '/ref/'
    if not os.path.exists(os.environ['CRDS_PATH']):
        utils.mkdir(os.environ['CRDS_PATH'])
        if glob('*flc.fits'):
            os.system('python -m crds.bestrefs --update-bestrefs --sync-references=1 --files *flc.fits')
        if glob('*flt.fits'):
            os.system('python -m crds.bestrefs --update-bestrefs --sync-references=1 --files *flt.fits')
        if glob('*c0m.fits'):
            os.system('python -m crds.bestrefs --update-bestrefs --sync-references=1 --files *c0m.fits')


def get_filter_name(image):
    if 'c0m.fits' in image:
        f = fits.getval(image, 'FILTNAM1')
        if len(f) == 0:
            f = fits.getval(image, 'FILTNAM2')
    else:
        try:
            f = fits.getval(image, 'FILTER')
        except KeyError:
            f = fits.getval(image, 'FILTER1')
            if 'clear' in f.lower():
                f = fits.getval(image, 'FILTER2')
    return f.lower()


def get_raw_image_filenames(raw_directory):
    images = []
    for ftype in ['flc', 'flt', 'c0m']:
        img_base = [img.split('_')[0] for img in images]
        ftype_images = glob(os.path.join(raw_directory, '*{0}.fits'.format(ftype)))
        for img in ftype_images:
            if img.split('_')[0] not in img_base:
                images.append(img)
    return images


def get_default_reference_hst(visit_meta_data):
    preferred_filters = ['f555w', 'f625w', 'f606w', 'f622w', 'f814w',
                         'f791w', 'f438w', 'f435w', 'f110w', 'f160w']
    for filt in preferred_filters:
        if filt in visit_meta_data['filter']:
            visit = visit_meta_data[visit_meta_data['filter'] == filt][0]
            break
    else:  # choose the central filter
        visit_meta_data.sort('filter')
        visit = visit_meta_data[len(visit_meta_data) // 2]
    return visit


def get_target_coordinates(flt_file):
    hdulist = fits.open(flt_file)
    ra = hdulist[0].header['RA_TARG']
    dec = hdulist[0].header['DEC_TARG']
    return [ra, dec]


def sort_raw_data(images, min_visit_separation=0.2):
    metadata = table.Table()
    metadata['filename'] = images
    metadata['expstart'] = [fits.getval(image, 'EXPSTART') for image in images]
    metadata['instrument'] = [fits_utils.get_instrument(image) for image in images]
    metadata['visit'] = ['visit1' for image in images]  # initialize the column
    metadata['filter'] = [get_filter_name(image) for image in images]
    metadata.sort('expstart')
    visit_starts, = np.where(np.diff(metadata['expstart']) > min_visit_separation)
    for visit_num, first_row in enumerate(visit_starts):
        metadata['visit'][first_row + 1:] = 'visit{:d}'.format(visit_num + 2)

    for image in metadata:
        folder_name = os.path.join(image['instrument'], image['visit'], image['filter'])
        utils.copy_if_not_exists(image['filename'], folder_name)

    visit_list = metadata.group_by(['instrument', 'visit', 'filter'])

    return visit_list.groups.keys


if __name__ == '__main__':
    run()
