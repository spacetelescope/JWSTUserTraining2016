import sl_accessor_funcs as funcs


def get_stray_light_bg(ra, dec, date, ra_dec_str, date_str):

    return funcs.get_sl(ra, dec, date, ra_dec_str, date_str)
