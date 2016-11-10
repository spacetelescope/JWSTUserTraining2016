'''
Script to create scamp dispersion reference files for NIRISS from input .csv resolving powers.
A constant pixel sampling is assumed.

'''
import numpy as np
import scipy.interpolate as ip
from astropy.io import fits
import astropy.io.ascii as at

# Actual sampling is currently not available. Assume Nyquist for now.
sampling = 2.
npix = 2048

files = [{'input': 'GR150C_R.csv', 'output': 'GR150C_disp.fits', 'sampling': 2.},
         {'input': 'GR150R_R.csv', 'output': 'GR150R_disp.fits', 'sampling': 2.},
         {'input': 'GR700XD_1_R.csv',
          'output': 'GR700XD_1_disp.fits',
          'sampling': 2.},
         {'input': 'GR700XD_2_R.csv', 'output': 'GR700XD_2_disp.fits', 'sampling': 2.}]

path = '../etc3d/refdata/NIRISS/dispersion/'

for file in files:
    sampling = file['sampling']
    data = at.read(path + file['input'])
    wave = data['col1']
    dlds = ip.interp1d(wave, wave / data['col2'] / sampling)
    rpower = ip.interp1d(wave, data['col2'])
    lrange = (np.min(wave), np.max(wave))

    wave_pt = lrange[0]
    wave_pix = []
    while wave_pt < lrange[1]:
        wave_pix.append(wave_pt)
        wave_pt += dlds(wave_pt)

    wave_pix = np.array(wave_pix)
    npix = wave_pix.shape[0]
    pixels = np.arange(npix)

    c1 = fits.Column(name='PIXELS', unit='PIXELS', format='1D', disp='I4', array=pixels)
    c2 = fits.Column(name='WAVELENGTH', unit='UM', format='1D', disp='F10.1', array=wave_pix)
    c3 = fits.Column(name='DLDS', unit='UM', format='1D', disp='G12.5', array=dlds(wave_pix))
    c4 = fits.Column(name='R', unit='TRANSMISSION', format='1D', disp='G12.5', array=rpower(wave_pix))

    coldefs = fits.ColDefs([c1, c2, c3, c4])

    hdu = fits.PrimaryHDU()
    tbhdu = fits.BinTableHDU.from_columns(coldefs)
    thdulist = fits.HDUList([hdu, tbhdu])

    thdulist.writeto(path + file['output'], clobber=True)
