import numpy as np
from astropy.io import fits
import astropy.io.ascii as at


def csv2fits_thruput(input, output, wave_scale=1., trans_scale=1.):

    data = at.read(input, comment='#')
    wave = data['col1']
    throughput = data['col2']

    c1 = fits.Column(name='WAVELENGTH', unit='UM', format='1D', disp='F10.1', array=wave)
    c2 = fits.Column(name='THROUGHPUT', unit='TRANSMISSION', format='1D', disp='G12.5', array=throughput)

    coldefs = fits.ColDefs([c1, c2])

    hdu = fits.PrimaryHDU()
    tbhdu = fits.BinTableHDU.from_columns(coldefs)
    thdulist = fits.HDUList([hdu, tbhdu])

    thdulist.writeto(output, clobber=True)
