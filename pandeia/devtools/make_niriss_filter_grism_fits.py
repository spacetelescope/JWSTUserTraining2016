#!/usr/bin/env python
#
import os, sys
import pysynphot as psyn

TOP = '/Users/sontag/dev/pandeia_data/jwst/niriss'
glist = ['gr150r', 'gr150c']
flist = ['f090w', 'f115w', 'f140m', 'f150w', 'f158m', 'f200w']

for grism in glist:
    for filter in flist:
        gname = TOP+'/blaze/'+grism+'_1_blaze.fits'
        gcurve = psyn.FileBandpass(gname)
        gcurve.convert('angstroms')

        fname = TOP+'/filters/'+filter+'.fits'
        fcurve = psyn.FileBandpass(fname)
        fcurve.convert('angstroms')

        outname = '/Users/sontag/NEW/niriss_'+grism+'_'+filter+'.fits'
        combined = fcurve*gcurve
        combined.convert('angstroms')
        print '-'*75
        print outname
        print 'fcurve (waveunits, wave.max, tp.max)'
        print fcurve.waveunits, fcurve.wave.max(), fcurve.throughput.max()
        print 'gcurve (waveunits, wave.max, tp.max)'
        print gcurve.waveunits, gcurve.wave.max(), gcurve.throughput.max()
        print 'combined (waveunits, wave.max, tp.max)'
        print combined.waveunits, combined.wave.max(), combined.throughput.max()
        print '-'*75
        combined.writefits(outname, precision='d')
