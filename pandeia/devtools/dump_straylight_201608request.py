#!/usr/bin/env python
#
import sys
from pandeia.engine.helpers.straylight_accessor import get_stray_light_bg


# Email from Jane R:
#
# Vicki,
# Yes, it would be extremely useful to have such a tarball, so that I can examine
# the contribution from each part of the background. Here are the fields of interest.
# 
# coords are celestial, J2000
# FIELD      RA            DEC             DATE
# HUDF       03 32 39.00   -27 47 29.0     Nov 1
# Lockman    10 45 00.00   58 00 00.0      Feb 8
# WC1.2zody  17 25 50.73   -73 19 35.99    please use the middle of visibility window
# 
# Do you need any other input  to proceed?
# Best,
# Jane R


# CONVERSIONS:
#
# HUDF:
#"ra_dec_str":"03 32 39.00  -27 47 29.0",
#"ra":53.1625,
#"dec":-27.79138888888889,
#"date":7244,
#"date_str":"Nov 1, 2019"
#
# Lockman:
#"ra_dec_str":"10 45 00.00   58 00 00.0"
#"ra":161.25,
#"dec":58,
#"date":6978,
#"date_str":"Feb 8, 2019"
#
# WC1.2zody:
#"ra_dec_str":"17 25 50.73   -73 19 35.99"
#"ra":261.461375
#"dec":-73.32666388888889
#"date":7109
#"date_str":"Jun 19, 2019"


def dump_sl_to_file(field_name, ra, dec, date, ra_dec_str, date_str):
    # get sl array
    sl_wave, sl_bg = get_stray_light_bg(ra, dec, date, ra_dec_str, date_str)
    assert len(sl_wave) == len(sl_bg), '?'
    # dump to file
    fname = field_name+'.dat'
    f = open(fname, 'w')
    for i in range(len(sl_wave)):
        f.write('%.4e  %.4e\n' % (sl_wave[i], sl_bg[i]))
    f.close()
    print('Dumped '+field_name+' to: '+fname)


#
# main routine
#
if __name__=='__main__': # in case something else imports this file

    dump_sl_to_file('HUDF',      53.1625,     -27.791389,  7244, '03 32 39.00  -27 47 29.00', 'Nov 1,  2019')
    dump_sl_to_file('Lockman',   161.25,      58,          6978, '10 45 00.00   58 00 00.00', 'Feb 8,  2019')
    dump_sl_to_file('WC1.2zody', 261.461375,  -73.3266639, 7109, '17 25 50.73  -73 19 35.99', 'Jun 19, 2019')

    sys.exit(0)
