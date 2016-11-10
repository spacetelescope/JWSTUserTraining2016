from shutil import copyfile

nircam_path = '../etc3d/refdata/NIRCam/filters/'
niriss_path = '../etc3d/refdata/NIRISS/filters/'

filters = ['F090W.fits', 'F115W.fits', 'F140M.fits', 'F150W.fits', 'F200W.fits', 'F277W.fits', 'F444W.fits',
           'F356W.fits', 'F430M.fits', 'F480M.fits']

for filter in filters:
    copyfile(nircam_path + filter, niriss_path + filter)
