#!/usr/bin/env python
import os, sys
import matplotlib
matplotlib.use("agg")
from matplotlib import pylab
import pyfits as pf

def create_plot(filename, title):
    """
    create_plot: Given a filename, plot the filter bandpass.

    Assume FITS table input.  Two columns. Column 1 is wavelength.
    Column 2 is throughput.

    Designed to work with the FORMAT of filter throughput files contained
    within SCAMP.


    arguments
    ---------
    filename: filename of FITS containing filter information
    title: String, representing the name desired for the plot

    """

    # Open the filter file
    hdul = pf.open(filename)

    # Extract the wavelength and throughput information
    wavelength = hdul[1].data['WAVELENGTH']
    throughput = hdul[1].data['THROUGHPUT']

    # Make plot
    pylab.plot(wavelength, throughput)
    pylab.xlabel(hdul[1].header['TUNIT1'])
#   pylab.xlabel("microns")
    pylab.ylabel("filter throughput")
    pylab.title(title)
    pylab.savefig(os.path.splitext(filename)[0]+'.png')
    pylab.close()

    #Close the open FITS file
    hdul.close()

# main
for fname in sys.argv[1:]:
    print os.path.splitext(fname)[0]+'.png'
    create_plot(fname, fname)
