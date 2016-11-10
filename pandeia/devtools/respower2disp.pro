PRO respower2disp, Rfile, dispfile, sampling=sampling

IF ~KEYWORD_SET(sampling) THEN sampling = 2.0

READCOL, Rfile, wave, R

dlds = wave/(R*sampling)

MWRFITS, dummy, dispfile+'.fits', /CREATE
MWRFITS, {wave:wave,dlds:dlds, waveunit:'micron'}, dispfile+'.fits'


END  
