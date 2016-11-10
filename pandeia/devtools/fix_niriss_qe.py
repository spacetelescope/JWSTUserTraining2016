import os
import astropy.io.fits as fits

path = '/Users/pontoppi/STSCI/SoftwareDevelopment/pandeia_data/jwst/detector'
files = [
         'jwst_niriss_h2rg_qe.fits'
         ]
         
for file in files:
    in_hdu = fits.open(os.path.join(path,file))
    wave = in_hdu[1].data['wavelength']
    throughput = in_hdu[1].data['qeff']
    
    w_col = fits.Column(name='WAVELENGTH', format='1D', unit='UM', disp="F10.1", array=wave)
    thru_col = fits.Column(name='THROUGHPUT', format='1D', unit='TRANSMISSION', disp="G12.5", array=throughput)

    cols = fits.ColDefs([w_col, thru_col])

    tb_hdu = fits.new_table(cols)
    hdr = fits.Header()

    pri_hdu = fits.PrimaryHDU(header=hdr)
    out_hdu = fits.HDUList([pri_hdu, tb_hdu])
    
    out_hdu.writeto(os.path.join(path,file),clobber=True)
