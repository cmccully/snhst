from astropy.io import fits
from astroscrappy import detect_cosmics
import logging
import numpy as np

from snhst.fits_utils import is_sci_extension

logger = logging.getLogger(__name__)


def detect_cosmic_rays(filename, options, output_image=None, masked_value=None):
    logger.info('Detecting cosmic rays in {f}'.format(f=filename))
    hdulist = fits.open(filename, mode='readonly')
    for hdu in hdulist:
        if is_sci_extension(hdu):
            if masked_value is not None:
                mask = hdu.data == masked_value
            else:
                mask = np.zeros(hdu.data.shape, dtype=np.bool)

            _crmask, crclean = detect_cosmics(hdu.data.copy().astype('<f4'), inmask=mask,
                                              readnoise=options['rdnoise'], gain=options['gain'],
                                              satlevel=options['saturation'],
                                              sigclip=options['sig_clip'],
                                              sigfrac=options['sig_frac'], objlim=options['obj_lim'])
            hdu.data[:, :] = crclean[:, :]

    if output_image is None:
        output_image = filename

    hdulist.writeto(output_image, overwrite=True)
    hdulist.close()
