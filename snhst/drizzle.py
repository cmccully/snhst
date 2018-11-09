from drizzlepac import astrodrizzle
from glob import glob
import numpy as np
from astropy.io import fits
from snhst import parameters
from snhst.find_offsets import run_tweakreg
from snhst.cosmic_rays import detect_cosmic_rays
from snhst import fits_utils
import os


def drizzle(images, output_name, options):

    output_name = os.path.join(os.path.dirname(images[0]), output_name)
    if '_flc.fits' in images[0]:
        output_filename = output_name + '_drc.fits'
    else:
        output_filename = output_name + '_drz.fits'

    if not os.path.exists(output_filename):
        # Based on the instrument, get the default parameter values
        # Keep any parameters that were provided by the user
        drizzle_options = parameters.get_instrument_parameters(images[0], options, 'drizzle')
        drizzle_options['output'] = output_name

        # If we want to find shifts, find the shifts comparing to the template image
        # (either the first image of the set or the provided template file)
        if drizzle_options['use_tweakshifts']:
            run_tweakreg(images, drizzle_options)
            drizzle_options['wcskey'] = 'TWEAK'

        # run astrodrizzle on the images
        run_astrodrizzle(images, drizzle_options)

        # Clean the cosmic rays from the image
        clean_cosmic_rays(output_filename, drizzle_options)

    return output_filename


def run_astrodrizzle(input_files, drizzle_options):

    input_dict = {'runfile': '', 'build': True, 'preserve': False,
                  'skystat': 'mode', 'skylower': 0.,
                  'driz_sep_fillval': -100000., 'driz_sep_wcs': True,
                  'combine_maskpt': 0.2, 'combine_lthresh': -10000.,
                  'combine_type': 'mean' if len(input_files) == 1 else 'minmed',
                  'final_fillval': -50000., 'final_units': 'counts', 'final_wcs': True}

    # These parameters were either set by default or passed in by the user
    for key in ['output', 'final_pixfrac', 'clean', 'num_cores', 'skysub']:
        if key in drizzle_options:
            input_dict[key] = drizzle_options[key]

    if 'refimage' in drizzle_options:
        driz_opts_to_include = ['bits', 'refimage']
    else:
        driz_opts_to_include = ['bits', 'rot', 'scale', 'outnx', 'outny', 'ra', 'dec']

    for key in driz_opts_to_include:
        if key in drizzle_options:
            input_dict['driz_sep_' + key] = drizzle_options[key]
            input_dict['final_' + key] = drizzle_options[key]

    astrodrizzle.AstroDrizzle(input_files, **input_dict)

    # add the sky value back in
    # only add the sky to non flagged pixels (-50k)
    if '_flc.fits' in input_files[0]:
        output_filename = drizzle_options['output'] + '_drc.fits'
    else:
        output_filename = drizzle_options['output'] + '_drz.fits'

    sci_hdu = fits.open(output_filename)
    sci_data = sci_hdu['SCI'].data
    mask = sci_data <= input_dict['final_fillval']
    if drizzle_options.get('skysub') in [True, None]:
        sci_data[~mask] -= np.min(sci_data[~mask])
    sci_data[mask] = 0.0

    sci_hdu.writeto(output_filename, overwrite=True)

    # remove temporary files missed by astrodrizzle's clean = True
    for filename in glob('*_skymatch_mask*.fits') + glob('*_staticMask.fits'):
        os.remove(filename)


def clean_cosmic_rays(image, options):
    if fits_utils.get_instrument(image) != 'wfc3_ir_full':
        output_image = image.replace('_dr', '_cr')
        detect_cosmic_rays(image, options['crpars'], output_image=output_image, masked_value=0.0)
