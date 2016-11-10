PRO miri_tframe, xstart, ystart, xsize, ysize

readmode = 1.
rowresets = 3.

colstart = (xstart-1)/4+1
colstop  = (xstart+xsize-2)/4+1
rowstart = ystart
rowstop  = ystart-1.+ysize

IF colstart EQ 1 THEN BEGIN 
   left_ref_pix = 3
ENDIF ELSE BEGIN
   left_ref_pix = 0
ENDELSE

IF colstop EQ 258 THEN BEGIN
   right_ref_pix = 3
ENDIF ELSE BEGIN
   right_ref_pix = 0
ENDELSE

roi_row_clks = (colstop-colstart+1.+left_ref_pix+right_ref_pix)*readmode
ovr_row_clks = colstart+rowresets+3.

row_clks = roi_row_clks+ovr_row_clks
rst_clks = 2.+rowresets

frame_clks = (rowstart-1.)*rst_clks + (1024.-rowstop)*rst_clks + (rowstop-rowstart+1.)*row_clks

frame_time = frame_clks*10e-6 ;seconds

print, 'frame clocks:', frame_clks
print, 'frame time:', frame_time, ' s'
END
