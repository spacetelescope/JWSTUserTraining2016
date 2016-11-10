#!/usr/bin/env python

import webbpsf as wp
import numpy as np
import os
import pyfits as pf


def psf_suite(nw=30, wmin=0.5, wmax=6., instrument=wp.NIRCam(),
              outname='PSF', fov_arcsec=1.5, aperture='any'):

    nw = float(nw)
    waves = (np.arange(nw) / (nw - 1)) ** 2. * (wmax - wmin) + wmin
    waves_m = waves * 1e-6

    for wave in waves_m:
        longname = outname + '_{0:.4f}'.format(wave * 1e6) + '.fits'
        longname = longname.lower()
        if os.path.isfile(longname):
            os.remove(longname)
        instrument.calcPSF(oversample=2, fov_arcsec=fov_arcsec, monochromatic=wave,
                           outfile=longname, clobber='true')

        # Add mode keyword to header
        hdulist = pf.open(longname, mode='update')
        prihdr = hdulist[0].header
        prihdr['APERTURE'] = (aperture, 'The observing aperture within the instrument FOV')
        hdulist.flush()
        hdulist.close()

doMIRI = False
doNIRSpec = False
doNIRCam = False
doNIRISS = True

if doMIRI:
    MIRI_FQPM = wp.MIRI()
    MIRI_FQPM.pupil_mask = 'MASKFQPM'
    MIRI_FQPM.image_mask = 'FQPM1065'
    psf_suite(nw=10, wmin=9.5, wmax=11.5, instrument=MIRI_FQPM,
              outname='MIRI_FQPM1065', fov_arcsec=12, aperture='FQPM1065')
    MIRI_FQPM.image_mask = 'FQPM1140'
    psf_suite(nw=10, wmin=10.5, wmax=12.5, instrument=MIRI_FQPM,
              outname='MIRI_FQPM1140', fov_arcsec=12, aperture='FQPM1140')
    MIRI_FQPM.image_mask = 'FQPM1550'
    psf_suite(nw=10, wmin=14.5, wmax=16.5, instrument=MIRI_FQPM,
              outname='MIRI_FQPM1550', fov_arcsec=12, aperture='FQPM1550')
    MIRI_FQPM.pupil_mask = 'MASKLYOT'
    MIRI_FQPM.image_mask = 'LYOT2300'
    psf_suite(nw=10, wmin=19., wmax=26.5, instrument=MIRI_FQPM,
              outname='MIRI_LYOT2300', fov_arcsec=12, aperture='LYOT2300')

    MIRI = wp.MIRI()
    MIRI.pixelscale = 0.196
    psf_suite(nw=30, wmin=4., wmax=8., instrument=MIRI,
              outname='MIRI_CH1', fov_arcsec=3.7, aperture='ch1')
    MIRI = wp.MIRI()
    MIRI.pixelscale = 0.196
    psf_suite(nw=30, wmin=7., wmax=12., instrument=MIRI,
              outname='MIRI_CH2', fov_arcsec=4.5, aperture='ch2')
    MIRI = wp.MIRI()
    MIRI.pixelscale = 0.245
    psf_suite(nw=30, wmin=11., wmax=19., instrument=MIRI,
              outname='MIRI_CH3', fov_arcsec=6.1, aperture='ch3')
    MIRI = wp.MIRI()
    MIRI.pixelscale = 0.273
    psf_suite(nw=30, wmin=17., wmax=30., instrument=MIRI,
              outname='MIRI_CH4', fov_arcsec=7.7, aperture='ch4')
    MIRI = wp.MIRI()
    MIRI.pixelscale = 0.11
    psf_suite(nw=30, wmin=4., wmax=30., instrument=MIRI,
              outname='MIRI_Imager', fov_arcsec=6, aperture='Imager')

    #MIRI = wp.MIRI()
    #MIRI.pixelscale = 0.11
    #MIRI.pupil_mask = 'P750L LRS grating'
    # psf_suite(nw=30,wmin=4.,wmax=15.,instrument=MIRI,outname='MIRI_LRS',fov_arcsec=5.5,aperture='LRSSLIT')

if doNIRCam:
    # Set these up to read in the current JWST.cfg to get the appropriate pixel sizes, rather than hard coding them.
    NIRCam = wp.NIRCam()
    NIRCam.pixelscale = 0.032
    psf_suite(nw=30, wmin=0.5, wmax=6.0, instrument=NIRCam,
              outname='NIRCam_SW', fov_arcsec=2., aperture='SW')

    NIRCam = wp.NIRCam()
    NIRCam.pixelscale = 0.064
    psf_suite(nw=30, wmin=0.5, wmax=6.0, instrument=NIRCam,
              outname='NIRCam_LW', fov_arcsec=2., aperture='LW')

if doNIRSpec:
    NIRSpec = wp.NIRSpec()
    NIRSpec.pixelscale = 0.105
    psf_suite(
        nw=30,
        wmin=0.5,
        wmax=6.0,
        instrument=NIRSpec,
        outname='NIRSpec_SLIT',
        fov_arcsec=3.0,
        aperture='Shutter,A200,A400,A1600')

    NIRSpec = wp.NIRSpec()
    NIRSpec.pixelscale = 0.105
    psf_suite(nw=30, wmin=0.5, wmax=6.0, instrument=NIRSpec,
              outname='NIRSpec_IFU', fov_arcsec=3.0, aperture='IFU')

if doNIRISS:
    NIRISS = wp.NIRISS()
    NIRISS.pixelscale = 0.0656
    psf_suite(nw=30, wmin=0.5, wmax=6.0, instrument=NIRISS, outname='NIRISS', fov_arcsec=2., aperture='Imager')
#    NIRISS = wp.NIRISS()
#    NIRISS.pupil_mask = 'GR700XD'
#    NIRISS.pixelscale = 0.0656
#    psf_suite(nw=20, wmin=0.6, wmax=2.6, instrument=NIRISS, outname='NIRISS_SOSS', fov_arcsec=2., aperture='GR700XD')
