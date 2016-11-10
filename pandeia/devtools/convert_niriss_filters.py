from csv2fits import csv2fits_thruput as csv2fits

infiles = [
    '../etc3d/refdata/NIRISS/filters/F158M.csv',
    '../etc3d/refdata/NIRISS/filters/F380M.csv']
outfiles = [
    '../etc3d/refdata/NIRISS/filters/F158M.fits',
    '../etc3d/refdata/NIRISS/filters/F380M.fits']

for infile, outfile in zip(infiles, outfiles):
    csv2fits(infile, outfile)
