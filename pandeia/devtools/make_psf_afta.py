#!/usr/bin/env python

from webbpsf import wfirst
import numpy as np
import os
from astropy.io import fits

def psf_suite(nw=30, wmin=0.5, wmax=6., instrument=wfirst.WFIRSTImager(),
              outname='PSF', fov_arcsec=4.0, aperture='any', alt_instrument=None):

    nw = float(nw)
    waves = (np.arange(nw) / (nw - 1)) ** 2. * (wmax - wmin) + wmin
    waves_m = waves * 1e-6

    for wave in waves_m:
        longname = outname + '_{0:.4f}'.format(wave * 1e6) + '.fits'
        if os.path.isfile(longname):
            os.remove(longname)
        instrument.calcPSF(oversample=4, fov_arcsec=fov_arcsec, monochromatic=wave,
                           outfile=longname, clobber='true')

        # Add mode keyword to header
        hdulist = fits.open(longname, mode='update')
        prihdr = hdulist[0].header
        prihdr['APERTURE'] = (aperture, 'The observing aperture within the instrument FOV')
        if alt_instrument is not None:
            prihdr['INSTRUME'] = (alt_instrument, '')
            
        
        hdulist.flush()
        hdulist.close()

doWFIRSTImager = False
doWFIRSTIFU = True

if doWFIRSTImager:
    WFIRSTImager = wfirst.WFIRSTImager()
    psf_suite(nw=30, wmin=0.4, wmax=2.6, instrument=WFIRSTImager,
              outname='WFIRSTImager_any', fov_arcsec=5)

if doWFIRSTIFU:
    WFIRSTImager = wfirst.WFIRSTImager()
    psf_suite(nw=30, wmin=0.4, wmax=2.6, instrument=WFIRSTImager,
              outname='WFIRSTIFU_ifu', fov_arcsec=5, aperture='ifu', alt_instrument='wfirstifu')
