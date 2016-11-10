PRO cdbs2scamp, in, out, wave_scale=wave_scale, trans_scale=trans_scale

IF ~KEYWORD_SET(wave_scale) THEN wave_scale = 1e4 ;angstrom->micron
IF ~KEYWORD_SET(trans_scale) THEN trans_scale = 1.

hdu_main = MRDFITS(in,0,HDR)
table    = MRDFITS(in,1)

MWRFITS, dummy, out, HDR,/CREATE
MWRFITS, {wave:table.wavelength/wave_scale,transmission:table.throughput/trans_scale}, out


END
