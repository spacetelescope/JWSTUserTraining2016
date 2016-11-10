PRO create_NIRSpec_dispersion, wavecal_dets, R_table, outfile

;wavecal_dets = filenames for the two SCAs. 
dlds  = []
wave_all = []
xsize = 2048
gap   = 148
xdet  = FINDGEN(xsize)
xgap  = FINDGEN(gap)

xall = [xdet,xgap+xsize,xdet+xgap+xsize]

FOR i=0,N_ELEMENTS(wavecal_dets)-1 DO BEGIN
   wavecal = MRDFITS(wavecal_dets[i],1)
                                ;There are slight differences between
                                ;the apertures (as there will be
                                ;across the MSA). We do not model that
                                ;at the moment, and are likely to
                                ;never include this complexity for NIRSpec.
   dlds    = [dlds,DERIV(wavecal.s,wavecal.a200_2)]
   wave_all = [wave_all,wavecal.a200_2]
ENDFOR

READCOL, R_table, wR, R

R_int = INTERPOL(R,wR,wave_all)

MWRFITS, dummy, outfile, /CREATE
MWRFITS, {pixels:xall,wave:wave_all,dlds:dlds,R:R_int, wave_unit:'microns'}, outfile, /CREATE


END


PATH    = '/Users/pontoppi/STSCI/SoftwareDevelopment/scamp/etc3d/refdata/NIRSpec/wavecal/'
OUTPATH = '/Users/pontoppi/STSCI/SoftwareDevelopment/scamp/etc3d/refdata/NIRSpec/dispersion/'
waveval_dets_G140H = PATH+['G140H_SCA491_1Dwave.fits','G140H_SCA492_1Dwave.fits']
waveval_dets_G235H = PATH+['G235H_SCA491_1Dwave.fits','G235H_SCA492_1Dwave.fits']
waveval_dets_G395H = PATH+['G395H_SCA491_1Dwave.fits','G395H_SCA492_1Dwave.fits']
waveval_dets_G140M = PATH+['G140M_SCA491_1Dwave.fits','G140M_SCA492_1Dwave.fits']
waveval_dets_G235M = PATH+['G235M_SCA491_1Dwave.fits','G235M_SCA492_1Dwave.fits']
waveval_dets_G395M = PATH+['G395M_SCA491_1Dwave.fits','G395M_SCA492_1Dwave.fits']

R_G140H = OUTPATH+'G140H_R.csv'
R_G235H = OUTPATH+'G235H_R.csv'
R_G395H = OUTPATH+'G395H_R.csv'
R_G140M = OUTPATH+'G140M_R.csv'
R_G235M = OUTPATH+'G235M_R.csv'
R_G395M = OUTPATH+'G395M_R.csv'


create_NIRSPEC_dispersion, waveval_dets_G140H, R_G140H, OUTPATH+'G140H_disp.fits'
create_NIRSPEC_dispersion, waveval_dets_G235H, R_G235H, OUTPATH+'G235H_disp.fits'
create_NIRSPEC_dispersion, waveval_dets_G395H, R_G395H, OUTPATH+'G395H_disp.fits'
create_NIRSPEC_dispersion, waveval_dets_G140M, R_G140M, OUTPATH+'G140M_disp.fits'
create_NIRSPEC_dispersion, waveval_dets_G235M, R_G235M, OUTPATH+'G235M_disp.fits'
create_NIRSPEC_dispersion, waveval_dets_G395M, R_G395M, OUTPATH+'G395M_disp.fits'

END
