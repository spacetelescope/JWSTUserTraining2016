import os
import astropy.io.fits as fits
import numpy as np

path = '/Users/pontoppi/STSCI/SoftwareDevelopment/pandeia_data/jwst/nirspec/blaze'
files = ['jwst_nirspec_g140h_speceff.fits',
         'jwst_nirspec_g140m_speceff.fits',
         'jwst_nirspec_g235h_speceff.fits',
         'jwst_nirspec_g235m_speceff.fits',
         'jwst_nirspec_g395h_speceff.fits',
         'jwst_nirspec_g395m_speceff.fits',
         'jwst_nirspec_mirror_speceff.fits',
         'jwst_nirspec_prism_speceff.fits']
         
for file in files:
    fullpath = os.path.join(path,file)
    fitsfile = fits.open(fullpath,mode='update')
    
    # Remove the 10% artificial margin the files were delivered with.
    fitsfile[1].data['throughput'] *= 1.1
    
    # Modify the header to match the content:
    fitsfile[0].header['DESCRIP'] = 'GWA element efficiency - as-built'
    fitsfile[0].header.add_history('K. M. Pontoppidan: Removed 10 percent margin; margin should not be hidden in as-built reference files.')
    
    fitsfile.flush()
    fitsfile.close()