PRO csv2fits_background, file, wave_scale=wave_scale, trans_scale=trans_scale,  ext=ext

IF ~KEYWORD_SET(wave_scale) THEN wave_scale = 1.
IF ~KEYWORD_SET(trans_scale) THEN trans_scale = 1.
IF ~KEYWORD_SET(ext) THEN ext = 'csv'

readcol, file+'.'+ext, x, y
MWRFITS, dummy, file+'.fits', /CREATE
MWRFITS, {wave:x/wave_scale,sb:y/trans_scale,sb_unit:'MJy/sr'}, file+'.fits'

END
