import pyfits as pf
import numpy as np

wave = np.linspace(0.3,1.1,100)
transmission = wave*0+0.9

c1  = pf.Column(name='wave', format='D', array=wave)
c2  = pf.Column(name='transmission', format='D', array=transmission)
coldefs = pf.ColDefs([c1,c2])

tbhdu   = pf.new_table(coldefs)
hdu     = pf.PrimaryHDU()
thdulist = pf.HDUList([hdu,tbhdu])

filename = 'ccd_qe.fits'
thdulist.writeto(filename,clobber=True)

