from csv2fits import *

infiles = ['../etc3d/refdata/NIRISS/optical/clear.csv',
           '../etc3d/refdata/NIRISS/optical/clearp.csv',
           '../etc3d/refdata/NIRISS/optical/pom_and_tmas.csv']

outfiles = ['../etc3d/refdata/NIRISS/optical/clear.fits',
            '../etc3d/refdata/NIRISS/optical/clearp.fits',
            '../etc3d/refdata/NIRISS/optical/pom_and_tmas.fits']

for infile, outfile in zip(infiles, outfiles):
    csv2fits(infile, outfile)
