#!/usr/bin/env python

# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
web2fits.py

This script takes NIRSpec dispersion data as presented at:

http://www.stsci.edu/jwst/instruments/nirspec/sensitivity/

and creates a FITS file suitable for use by the JWST ETC engine.
"""
import sys
import os
import numpy as np
import astropy.io.fits as fits
import scipy.interpolate as ip

try:
    file = sys.argv[1]
    wave, r = np.loadtxt(file, unpack=True)
except Exception as e:
    print str(e)
    print "\nUsage: web2fits.py <columnar ascii data file>"
    sys.exit()

# need to get wave_pix from an existing file that's known to work
hdu = fits.open("G395H_disp.fits")
t = hdu[1]
wave_pix = t.data['WAVELENGTH']

# wave_pix extends beyond the range given in the data found on the web pages.  fortunately,
# the resolving power is a linear function of wavelength.
new_wave = newwave = np.arange(400, 1100)/200.0  # same sampling, but wider range
fit = np.polyfit(wave, r, deg=1)  # simple linear fit
new_r = fit[1] + fit[0] * new_wave

r_interp = ip.interp1d(new_wave, new_r)  # new interp1d based on new_wave, new_r
r_pix = r_interp(wave_pix)
dlds = 0.5 * wave_pix / r_pix

w_col = fits.Column(name='WAVELENGTH', format='1D', unit='UM', disp="F10.1", array=wave_pix)
disp_col = fits.Column(name='DLDS', format='1D', unit='UM', disp="G12.5", array=dlds)
r_col = fits.Column(name='R', format='1D', unit='TRANSMISSION', disp="G12.5", array=r_pix)

cols = fits.ColDefs([w_col, disp_col, r_col])

tb_hdu = fits.BinTableHDU.from_columns(cols)
hdr = fits.Header()

pri_hdu = fits.PrimaryHDU(header=hdr)
out_hdu = fits.HDUList([pri_hdu, tb_hdu])

outfile = "disp.fits"
out_hdu.writeto(outfile)
