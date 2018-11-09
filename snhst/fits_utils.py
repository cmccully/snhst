from astropy.io import fits


def is_sci_extension(hdu):
    return 'EXTNAME' in hdu.header and hdu.header['EXTNAME'].upper() == 'SCI'


def get_instrument(filename):
    hdu = fits.open(filename, mode='readonly')
    instrument = hdu[0].header['INSTRUME'].lower()
    if instrument.upper == 'WFPC2':
        detector = 'wfpc2'
        subarray = 'full'
    else:
        detector = hdu[0].header['DETECTOR'].lower()
        if hdu[0].header['SUBARRAY']:
            subarray = "sub"
        else:
            subarray = "full"
    return "{instrument}_{detector}_{subarray}".format(instrument=instrument, detector=detector,
                                                       subarray=subarray)


def get_data_extension(filename):
    """Get the science hdu of an image (may be compressed)"""
    reference_hdulist = fits.open(filename)
    for reference_hdu in reference_hdulist:
        if reference_hdu.header.get('NAXIS'):
            return reference_hdu
