#!/usr/bin/env python
#
import sys
from pandeia.engine.custom_exceptions import StraylightPositionError
from pandeia.engine.helpers.straylight_accessor import get_stray_light_bg


# Chris,
# Would it be possible to please run one more stray light case?
# The target is Sag A*, at 17h45m40.0s -29d00m28s.  Midpoint of the visibility is fine.
# Thanks!
# Jane R.


# CONVERSIONS:
#
# Sag A:
# radec = 17:45:40.0 -29:00:28.0 on Jan 1, 2019 ...
RA_DEC_STR = "17:45:40.0 -29:00:28.0"
RA = 266.4166666666667
DEC = -29.007777777778
#"date":6940
#"date_str":"Jan 1, 2019"

# for date Dec 31, 2019,
# "date":7304



def is_visible(field_name, ra, dec, date, ra_dec_str, date_str):
    # get sl array
    vis_state = "visible"
    try:
        sl_wave, sl_bg = get_stray_light_bg(ra, dec, date, ra_dec_str, date_str)
    except StraylightPositionError:
        vis_state = "NOT visible"
    if vis_state == "visible":
        assert len(sl_wave) == len(sl_bg), '?'
    return '%s is %s on MJD2K=%s (radec %s)' % (field_name, vis_state, date_str, ra_dec_str)


#
# main routine
#
if __name__=='__main__': # in case something else imports this file

    for mjd in range(6940, 7304+1):
        print is_visible("Sag A", RA, DEC, mjd, RA_DEC_STR, str(mjd))

    sys.exit(0)
