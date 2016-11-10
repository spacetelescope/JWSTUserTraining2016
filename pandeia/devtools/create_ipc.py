'''
Script to generate a 3x3 IPC kernel relevant for the NIR HAWAII 2RG JWST detectors.
'''

import astropy.io.fits as pf
import numpy as np

#Following the description of McCollough et al. 2007 (JWST-STScI-001053)
alpha = 0.026
beta = 0.015
ipc_kernel = np.array([
    [0,     beta,             0],
    [alpha, 1-2*alpha-2*beta, alpha],
    [0,     beta,             0]
])

hdu = pf.PrimaryHDU(ipc_kernel)
hdulist = pf.HDUList([hdu])

filename = 'HAWAII_2RG_IPC_kernel.fits'
hdulist.writeto(filename, clobber=True)
