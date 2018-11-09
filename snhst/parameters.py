from snhst import fits_utils

global_defaults = {
    'ground_reference': None,
    'ground_reference_catalog': None,
    'ra': None,
    'dec': None,
    'use_tweakshifts': True,
    'tweakshifts_threshold': 10,
    'use_sep': False,
    'make_visit_templates': True,
    'dolphot': {
        'FitSky': 2,
        'SkipSky': 2,
        'RCombine': 1.5,
        'SkySig': 2.25,
        'SecondPass': 5,
        'SigFindMult': 0.85,
        'MaxIT': 25,
        'NoiseMult': 0.10,
        'FSat': 0.999,
        'ApCor': 1,
        'RCentroid': 2,
        'PosStep': 0.25,
        'dPosMax': 2.5,
        'SigPSF': 10.0,
        'PSFres': 1,
        'Align': 4,
        'Rotate': 1,
        'ACSuseCTE': 0,
        'WFC3useCTE': 0,
        'WFPC2useCTE': 1,
        'FlagMask': 7,
        'SigFind': 2.5,
        'SigFinal': 3.5,
        'UseWCS': 1,
        'AlignIter': 3,
        'AlignTol': 0.5,
        'AlignStep': 0.2,
        'VerboseData': 1,
        'NegSky': 0,
        'Force1': 0,
        'DiagPlotType': 'PNG',
        'InterpPSFlib': 1,
        },
    'sep': {
        'thresh': 20.,
        'deblend_cont': 0.005,
        'minarea': 50.,
        },
    'fakelist': {
        'filter1': None,
        'filter2': None,  # None = set based on data filters
        'filter1_min': 20.,
        'filter1_max': 30.,
        'color_min': 0.,
        'color_max': 0.1,
        'nstar': 50000,
        },
    }

instrument_defaults = {
    'wfc3_uvis_full': {
        'input_files': '*_flc.fits',
        'drizzle': {
            'rot': 0.,
            'clean': True,
            'crpars': {
                'rdnoise': 6.5,
                'gain': 1.0,
                'saturation': 70000.0,
                'sig_clip': 4.0,
                'sig_frac': 0.2,
                'obj_lim': 6.0,
                },
            },
        'calcsky': {
            'r_in': 15,
            'r_out': 35,
            'step': 4,
            'sigma_low': 2.25,
            'sigma_high': 2.00,
            },
        'dolphot_img': {
            'apsky': '15 25',
            'RAper': 2,
            'RChi': 1.5,
            'RPSF': 10,
            'RSky': '15 35',
            'RSky2': '3 6',
            },
        },
    'wfc3_uvis_sub': {
        'input_files': '*_flt.fits',
        'drizzle': {
            'rot': 0.,
            'clean': True,
            'crpars': {
                'rdnoise': 6.5,
                'gain': 1.0,
                'saturation': 70000.0,
                'sig_clip': 4.0,
                'sig_frac': 0.2,
                'obj_lim': 6.0,
                },
            },
        'calcsky': {
            'r_in': 15,
            'r_out': 35,
            'step': 4,
            'sigma_low': 2.25,
            'sigma_high': 2.00,
            },
        'dolphot_img': {
            'apsky': '15 25',
            'RAper': 2,
            'RChi': 1.5,
            'RPSF': 10,
            'RSky': '15 35',
            'RSky2': '3 6',
            },
        },
    # Note that wfc3_ir does not ever do cosmic ray removal
    # because it automatically flags cosmic rays using up the ramp sampling
    'wfc3_ir_full': {
        'input_files': '*_flt.fits',
        'drizzle': {
            'rot': 0.,
            'clean': True,
            'crpars': {
                'rdnoise': 6.5,
                'gain': 1.0,
                'saturation': 70000.0,
                'sig_clip': 4.0,
                'sig_frac': 0.2,
                'obj_lim': 6.0,
                },
            },
        'calcsky': {
            'r_in': 10,
            'r_out': 25,
            'step': 2,
            'sigma_low': 2.25,
            'sigma_high': 2.00,
            },
        'dolphot_img': {
            'apsky': '15 25',
            'RAper': 2,
            'RChi': 1.5,
            'RPSF': 10,
            'SkipSky': 1,
            'RSky': '8 20',
            'RSky2': '3 6',
            },
        },
    'acs_wcs_full': {
        'input_files': '*_flc.fits',
        'drizzle': {
            'rot': 0.,
            'clean': True,
            'crpars': {
                'rdnoise': 6.5,
                'gain': 1.0,
                'saturation': 70000.0,
                'sig_clip': 3.0,
                'sig_frac': 0.1,
                'obj_lim': 5.0,
                },
            },
        'calcsky': {
            'r_in': 15,
            'r_out': 35,
            'step': 4,
            'sigma_low': 2.25,
            'sigma_high': 2.00,
            },
        'dolphot_img': {
            'apsky': '15 25',
            'RAper': 2,
            'RChi': 1.5,
            'RPSF': 10,
            'RSky': '15 35',
            'RSky2': '3 6',
            },
        },
    'acs_wcs_sub': {
        'input_files': '*_flt.fits',
        'drizzle': {
            'rot': 0.,
            'clean': True,
            'crpars': {
                'rdnoise': 6.5,
                'gain': 1.0,
                'saturation': 70000.0,
                'sig_clip': 3.0,
                'sig_frac': 0.1,
                'obj_lim': 5.0,
                },
            },
        'calcsky': {
            'r_in': 15,
            'r_out': 35,
            'step': 4,
            'sigma_low': 2.25,
            'sigma_high': 2.00,
            },
        'dolphot_img': {
            'apsky': '15 25',
            'RAper': 2,
            'RChi': 1.5,
            'RPSF': 10,
            'RSky': '15 35',
            'RSky2': '3 6',
            },
        },
    'wfpc2_wfpc2_full': {
        'input_files': '*_c0m.fits',
        'drizzle': {
            'rot': 0.,
            'clean': True,
            'crpars': {
                'rdnoise': 10.0,
                'gain': 7.0,
                'saturation': 27000.0,
                'sig_clip': 4.0,
                'sig_frac': 0.3,
                'obj_lim': 6.0,
                },
            },
        'calcsky': {
            'r_in': 10,
            'r_out': 25,
            'step': 2,
            'sigma_low': 2.25,
            'sigma_high': 2.00,
            },
        'dolphot_img': {
            'apsky': '15 25',
            'RAper': 2,
            'RChi': 1.5,
            'RPSF': 10,
            'RSky': '15 35',
            'RSky2': '3 6',
            },
        },
    }


def set_default_parameters(options, defaults):
    for key in defaults:
        if key not in options:
            options[key] = defaults[key]


def get_instrument_parameters(image, options, key):
    instrument_string = fits_utils.get_instrument(image)
    set_default_parameters(options, instrument_defaults[instrument_string][key])
    return options
