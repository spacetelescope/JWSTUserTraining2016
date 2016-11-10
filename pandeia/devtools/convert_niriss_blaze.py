from csv2fits import csv2fits_thruput as csv2fits

infiles = [
    '../etc3d/refdata/NIRISS/blaze/GR150C_0_blaze.txt',
    '../etc3d/refdata/NIRISS/blaze/GR150R_0_blaze.txt',
    '../etc3d/refdata/NIRISS/blaze/GR150C_1_blaze.txt',
    '../etc3d/refdata/NIRISS/blaze/GR150R_1_blaze.txt',
    '../etc3d/refdata/NIRISS/blaze/GR150C_2_blaze.txt',
    '../etc3d/refdata/NIRISS/blaze/GR150R_2_blaze.txt',
    '../etc3d/refdata/NIRISS/blaze/GR700XD_1_blaze.txt',
    '../etc3d/refdata/NIRISS/blaze/GR700XD_2_blaze.txt']
outfiles = [
    '../etc3d/refdata/NIRISS/blaze/GR150C_0_blaze.fits',
    '../etc3d/refdata/NIRISS/blaze/GR150R_0_blaze.fits',
    '../etc3d/refdata/NIRISS/blaze/GR150C_1_blaze.fits',
    '../etc3d/refdata/NIRISS/blaze/GR150R_1_blaze.fits',
    '../etc3d/refdata/NIRISS/blaze/GR150C_2_blaze.fits',
    '../etc3d/refdata/NIRISS/blaze/GR150R_2_blaze.fits',
    '../etc3d/refdata/NIRISS/blaze/GR700XD_1_blaze.fits',
    '../etc3d/refdata/NIRISS/blaze/GR700XD_2_blaze.fits']

for infile, outfile in zip(infiles, outfiles):
    csv2fits(infile, outfile)
