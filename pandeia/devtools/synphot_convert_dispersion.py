#!/usr/bin/env python

# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
import os

import astropy.io.fits as fits

try:
    file = sys.argv[1]
    in_hdu = fits.open(file)
except Exception as e:
    print str(e)
    print "\nUsage: synphot_convert.py <fitsfile>"
    sys.exit()

wave = in_hdu[1].data[0][1]
dlds = in_hdu[1].data[0][2]
r = in_hdu[1].data[0][3]

w_col = fits.Column(name='WAVELENGTH', format='1D', unit='UM', disp="F10.1", array=wave)
disp_col = fits.Column(name='DLDS', format='1D', unit='UM', disp="G12.5", array=dlds)
r_col = fits.Column(name='R', format='1D', unit='TRANSMISSION', disp="G12.5", array=r)

cols = fits.ColDefs([w_col, disp_col, r_col])

tb_hdu = fits.new_table(cols)
hdr = fits.Header()

pri_hdu = fits.PrimaryHDU(header=hdr)
out_hdu = fits.HDUList([pri_hdu, tb_hdu])

outfile = "new.fits"
out_hdu.writeto(outfile)
os.rename(outfile, file)
