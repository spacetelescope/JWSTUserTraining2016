import os
import astropy.io.fits as fits

path = '/Users/pontoppi/STSCI/SoftwareDevelopment/pandeia_data/jwst/niriss/dispersion'
files = ['jwst_niriss_gr150c-ordp1_disp.fits',
         'jwst_niriss_gr150c-ordp2_disp.fits',
         'jwst_niriss_gr150c-ordp3_disp.fits',
         'jwst_niriss_gr150r-ordp1_disp.fits',
         'jwst_niriss_gr150r-ordp2_disp.fits',
         'jwst_niriss_gr150r-ordp3_disp.fits',
         'jwst_niriss_gr700xd-ord1_disp.fits',
         'jwst_niriss_gr700xd-ord1_disp.fits',
         'jwst_niriss_gr700xd-ord2_disp.fits'
         ]
         
for file in files:
    in_hdu = fits.open(os.path.join(path,file))
    wave = in_hdu[1].data['wavelength']
    dlds = in_hdu[1].data['dispersion']
    r = in_hdu[1].data['resolv_power']
    
    w_col = fits.Column(name='WAVELENGTH', format='1D', unit='UM', disp="F10.1", array=wave)
    disp_col = fits.Column(name='DLDS', format='1D', unit='UM', disp="G12.5", array=dlds)
    r_col = fits.Column(name='R', format='1D', unit='TRANSMISSION', disp="G12.5", array=r)

    cols = fits.ColDefs([w_col, disp_col, r_col])

    tb_hdu = fits.new_table(cols)
    hdr = fits.Header()

    pri_hdu = fits.PrimaryHDU(header=hdr)
    out_hdu = fits.HDUList([pri_hdu, tb_hdu])
    
    out_hdu.writeto(os.path.join(path,file),clobber=True)
