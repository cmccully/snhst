import os

import numpy as np

from snhst import fits_utils, utils
from glob import glob
from astropy.io import fits
from astropy import table
from snhst import parameters
from reproject import reproject_interp
from snhst import reduce_hst_data
from snhst import wcs
from astropy.wcs import WCS
import re
import stwcs


def dolphot(template_image, images, path, options):
    # Copy the raw input files into a new directory
    images = [utils.copy_if_not_exists(image, path) for image in images]

    template_image = utils.copy_if_not_exists(template_image, path)

    for image in images:
        if needs_to_be_masked(image):
            mask_image(image)

        if needs_to_split_groups(image):
            split_groups(image)

    if needs_to_split_groups(template_image):
        split_groups(template_image)
    template_image = template_image.replace('.fits', '.chip1.fits')

    overlapping_images = get_overlapping_split_images(template_image, images)

    for image in overlapping_images:
        if needs_to_calc_sky(image):
            calc_sky(image, options.get('calcsky', {}))

    dp_out = os.path.join(path, 'dp.out')
    if not os.path.exists(dp_out):
        run_dolphot(template_image, overlapping_images, path, options)
    output_table = parse_dolphot_table(dp_out, os.path.join(path, 'hst.cat'))

    hdr = fits.getheader(template_image)
    output_table['ra'], output_table['dec'] = WCS(hdr).all_pix2world(output_table['x'], output_table['y'], 1)

    return output_table


def needs_to_be_masked(image):
    # Masking should remove all of the DQ arrays etc, so make sure that any extensions with data in
    # in them are only SCI extensions. This might not be 100% robust, but should be good enough.
    hdulist = fits.open(image)
    needs_masked = False
    for hdu in hdulist:
        if hdu.data is not None and 'EXTNAME' in hdu.header:
            if hdu.header['EXTNAME'].upper() != 'SCI':
                needs_masked = True
    return needs_masked


def mask_image(image):
    instrument = fits_utils.get_instrument(image).split('_')[0]
    os.system('{instrument}mask {image}'.format(instrument=instrument, image=image))


def needs_to_calc_sky(image):
    return not os.path.exists(image.replace('.fits', '.sky.fits'))


def calc_sky(image, options):
    calcsky_opts = parameters.get_instrument_parameters(image, options, 'calcsky')
    cmd = 'calcsky {image} {rin} {rout} {step} {sigma_low} {sigma_high}'
    cmd = cmd.format(image=image.replace('.fits', ''), rin=calcsky_opts['r_in'],
                     rout=calcsky_opts['r_out'], step=calcsky_opts['step'],
                     sigma_low=calcsky_opts['sigma_low'],
                     sigma_high=calcsky_opts['sigma_high'])
    print(cmd)
    os.system(cmd)


def needs_to_split_groups(image):
    return len(glob(image.replace('.fits', '.chip?.fits'))) == 0


def split_groups(image):
    os.system('splitgroups {filename}'.format(filename=image))


def get_overlapping_split_images(template_image, images):
    overlapping_images = []
    template_hdu = fits.open(template_image)
    for image in images:
        split_images = glob(image.replace('.fits', '.chip?.fits'))
        extensions_per_chip = {1: 1, 2: 4}
        for split_image in split_images:
            # strip off the .fits to get the chip number
            chip = int(split_image[-6])
            # Because of HST's new wcs convention we have to open a different file to get the wcs
            image_hdu = fits.open(image.replace('dolphot/', ''))
            image_wcs = stwcs.wcsutil.HSTWCS(image_hdu, extensions_per_chip[chip])
            # remap the template_image onto the image
            _, footprint = reproject_interp(template_hdu, image_wcs)
            # If the overlap is at least 10%, consider the image to be overlapping
            if (footprint > 0).sum() >= (0.1 * footprint.size):
                overlapping_images.append(split_image)
    return overlapping_images


def write_dolphot_image_parameters(file_object, image, i, options):
    file_object.write('img{i}_file = {file}\n'.format(i=i + 1, file=os.path.splitext(image)[0]))
    for par, value in parameters.get_instrument_parameters(image, options, 'dolphot_img').items():
        file_object.write('img{i}_{option} = {value}\n'.format(i=i + 1, option=par, value=value))


def write_dolphot_master_parameters(file_object, options):
    for par, value in options.items():
        file_object.write('{par} = {value}\n'.format(par=par, value=value))


def run_dolphot(template_image, images, path, options):
    f = open(os.path.join(path, 'dp.params'), 'w')
    parameters.set_default_parameters(options['dolphot'], parameters.global_defaults['dolphot'])
    write_dolphot_master_parameters(f, options['dolphot'])
    f.write('Nimg = {n}\n'.format(n=len(images)))
    f.write('img0_file = {drzfile}\n'.format(drzfile=os.path.splitext(template_image)[0]))
    for i, image in enumerate(images):
        write_dolphot_image_parameters(f, image, i, options.get('dolphot_img', {}))
    f.close()
    print(f'Running dolphot in {path}.')
    cmd = 'dolphot {0}/dp.out -p{0}/dp.params 2>&1 | tee -a {0}/dp.log'.format(path)
    print(f'{cmd}')
    dolphot_status = os.system(cmd)
    if dolphot_status:
        raise RuntimeError('Dolphot failed with status {}'.format(dolphot_status))


def cut_bad_dolphot_sources(catalog):
    # Reject bad sources
    catalog = catalog[(catalog['col11'] == 1) | (catalog['col11'] == 2)]
    # Sharpness cut
    catalog = catalog[np.abs(catalog['col7']) < 0.3]
    # Crowding cut
    catalog = catalog[catalog['col10'] < 0.5]
    return catalog


def parse_dolphot_table(hst_table, output_catalog):
    dolphot_catalog = table.Table.read(hst_table, format='ascii.fast_no_header')

    dolphot_catalog = cut_bad_dolphot_sources(dolphot_catalog)
    t = table.Table()
    t['x'] = dolphot_catalog['col3'] + 0.5
    t['y'] = dolphot_catalog['col4'] + 0.5
    t['mag'] = dolphot_catalog['col16']
    t['magerr'] = dolphot_catalog['col18']

    t.write(output_catalog, format='ascii.fast_no_header', overwrite=True)
    return t


def add_fake_stars(visit_meta_data, path, options):
    # find bluest and reddest filters
    filt_wl = np.array([re.search('[0-9]+', filt).group() for filt in visit_meta_data['filter']], float)
    filt_wl[filt_wl < 200.] *= 10.  # WFC3_IR filters are in different units
    if options['filter1'] is None:
        bluest_row = visit_meta_data[np.argmin(filt_wl)]
        inst = bluest_row['instrument'].split('_')[0]
        filt = bluest_row['filter']
        options['filter1'] = (inst + '_' + filt).upper()
    if options['filter2'] is None:
        reddest_row = visit_meta_data[np.argmax(filt_wl)]
        inst = reddest_row['instrument'].split('_')[0]
        filt = reddest_row['filter']
        options['filter2'] = (inst + '_' + filt).upper()

    # make list of fake stars for dolphot
    os.system('fakelist {path}/dp.out {filter1} {filter2} {filter1_min} {filter1_max}'
              ' {color_min} {color_max} -nstar={nstar} 2>&1 |'
              ' tee -a {path}/fakelist.out'.format(path=path, **options))

    # re-run dolphot with params file modified to include fake star list
    with open(os.path.join(path, 'dp.params'), 'a') as f:
        f.write('FakeStars = {}/fakelist.out\n'.format(path))

    os.system('dolphot {0}/dp.out -p{0}/dp.params 2>&1 | tee -a {0}/dp.log.fake'.format(path))
    output_cat = os.path.join(path, 'hst_fake.cat')
    parse_dolphot_table(os.path.join(path, 'dp.out.fake'), output_cat)
    return output_cat
