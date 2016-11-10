#!/usr/bin/env python

'''
This code was copied from the straylight area of the BMG svn repo (rev 352)

One minor modification was made:
    1) convert ra and dec from radians
'''

import numpy as np
import os

from pandeia.engine.custom_exceptions import StraylightPositionError
from pandeia.engine.custom_exceptions import StraylightDataError
from pandeia.engine.custom_exceptions import DatelessBGError

import accessor_globals as glb

from mod_healpix_func import ang2pix_ring

SL_NWAVE = len(glb.wavelist)
THERMAL_FNAME = 'thermal_curve_jwst_jrigby_1.0.csv'


def return_index(ra, dec):  # This takes radians!
    return ang2pix_ring(glb.nside, ra, dec)


nonzodi_pix_dtype = np.dtype(
    [
        ('pix_ra', 'f8'),
        ('pix_dec', 'f8'),
        ('upos', [('x', 'f8'), ('y', 'f8'), ('z', 'f8')]),
        ('nonzodi_bg', ('f8', SL_NWAVE)),
        ('iday_index', ('i4', glb.NUM_DAYS))
    ]
)

zodi_sl_dtype = np.dtype(
    [
        ('zodi_bg', ('f8', SL_NWAVE)),
        ('stray_light_bg', ('f8', SL_NWAVE))
    ]
)

dateless_index_dtype = np.dtype([('dateless_index', ('i4', 3))])


def get_sl(ra, dec, mjd2000, ra_dec_str, date_str):
    """
    Provides equivalent straylight background spectra for pointing and mjd2000.

    Parameters
    ----------
    ra : double
        Right Ascension [degrees]
    dec : double
        declination [degrees]
    mjd2000 : int
        Modified Julian Day 2000

    Returns
    -------
    tuple of:
    wave : numpy.ndarray
        Wavelengths of background values [microns]
    stray_light_bg : numpy.ndarray
        Equivalent in field background from the scattered zodi, ism, cbi, stellar.
    If the target is not in the Field of Regard on iday, then an exception is raised.
    """

    # Cache is based on DOY 2020 so convert MJD2000 to DOY.
    # 51544 is Jan 1, 2000
    iday = int(float(mjd2000) % 365.25) + 1

    # Find the HEALPix number and the subdirectory
    ipix = ang2pix_ring(glb.nside, ra*glb.D2R, dec*glb.D2R)
    ipix_dir = ipix // 100

    # Set up the numpy arrays
    wave = np.array(glb.wavelist, dtype='double')

    sl_path = os.environ['SL_CACHE_DIR']
    file_name = "%s/%04d/sl_pix_%06d.bin" % (sl_path, ipix_dir, ipix)
    if not os.path.exists(file_name):
        msg = 'Stray light data is not available for position (%s) on %s.' % (ra_dec_str, date_str)
        raise StraylightDataError(msg)
    nonzodi_bg = np.fromfile(file_name, dtype=nonzodi_pix_dtype, count=1)
    iday_pt = nonzodi_bg['iday_index'][0][iday-1]
    if iday_pt == -1:
        msg = 'Specified position (%s) is not observable on %s.  See the <a href="/doc/background_help.txt" target="_blank">docs</a> for more information.' % (ra_dec_str, date_str)
        raise StraylightPositionError(msg)
    zodi_sl_bgs = np.memmap(
        file_name,
        offset=nonzodi_pix_dtype.itemsize + zodi_sl_dtype.itemsize * iday_pt,
        mode='r',
        dtype=zodi_sl_dtype
    )
    stray_light_bg = np.array(zodi_sl_bgs['stray_light_bg'][0], dtype='double')
    return (wave, stray_light_bg)


def get_dateless_bg(ra, dec, level):
    """
    Provides equivalent straylight background spectra for pointing and day.

    Parameters
    ----------
    ra : double
        Right Ascension [degrees]
    dec : double
        declination [degrees]
    level : string
        Single character L,M,H

    Returns
    -------
    tuple of:
    doy : integer
        Day of year in 2020 where the backgrounds match the specified level.
        L,M,H = 10%, 50%, 90%.
    wave : numpy.ndarray
        Standard wavelengths of background values [microns]
    infield_bg : numpy.ndarray
        Infield background spectra from zodi, cib, and ism on Day of Year 2020 [MJy/str].
    stray_light_bg : numpy.ndarray
        Equivalent infield background sptectrum from stray light on Day of Year 2020 [MJy/str].
    thermal_wave : numpy.ndarray
        Wavelengths of equivilant thermal background values [microns].
    thermal_bg : numpy.ndarray
        Equivalent infield thermal background spectrum from JWST thermal self emission [MJy/str].
    """

    if level == 'L':
        ilevel = 0
    elif level == 'M':
        ilevel = 1
    elif level == 'H':
        ilevel = 2
    else:
        msg = 'Input level parameter, %s, is not L, M, or H.' % (level)
        raise DatelessBGError(msg)

    # Find the HEALPix number and the subdirectory
    ipix = ang2pix_ring(glb.nside, ra*glb.D2R, dec*glb.D2R)
    ipix_dir = ipix // 100

    # Set up the numpy arrays
    wave = np.array(glb.wavelist, dtype='double')

    sl_path = os.environ['SL_CACHE_DIR']
    file_name = "%s/dateless_background_index.bin" % (sl_path)
    if not os.path.exists(file_name):
        msg = 'dateless_background_index.bin file does not exist in %s.' % (sl_path)
        raise StraylightDataError(msg)
    dateless_index = np.memmap(
        file_name,
        offset=dateless_index_dtype.itemsize * ipix,
        mode='r',
        dtype=dateless_index_dtype
    )
    iday = dateless_index['dateless_index'][0][ilevel]
    file_name = "%s/%04d/sl_pix_%06d.bin" % (sl_path, ipix_dir, ipix)
    if not os.path.exists(file_name):
        msg = 'Dateless background data is not available for position (%s, %s).' % (ra, dec)
        raise DatelessBGError(msg)
    nonzodi_bg = np.fromfile(file_name, dtype=nonzodi_pix_dtype, count=1)
    iday_pt = nonzodi_bg['iday_index'][0][iday-1]
    # Get the spectra for the day
    zodi_sl_bgs = np.memmap(
        file_name,
        offset=nonzodi_pix_dtype.itemsize + zodi_sl_dtype.itemsize * iday_pt,
        mode='r',
        dtype=zodi_sl_dtype
    )
    nonzodi_bg = np.array(nonzodi_bg['nonzodi_bg'][0], dtype='double')
    zodi_bg = np.array(zodi_sl_bgs['zodi_bg'][0], dtype='double')
    infield_bg = nonzodi_bg + zodi_bg
    stray_light_bg = np.array(zodi_sl_bgs['stray_light_bg'][0], dtype='double')
    [thermal_wave, thermal_bg] = get_thermal_background()
    return (iday, wave, infield_bg, stray_light_bg, thermal_wave, thermal_bg)


def get_thermal_background():
    """
    Provides equivalent thermal background spectra. (renamed from get_thermal_bg to avoid confusion)

    Parameters
    ----------
    None

    Returns
    -------
    tuple of:
    wave : numpy.ndarray
        Wavelengths of background values [microns]
    thermal_bg : numpy.ndarray
        Equivalent in field background from observatory thermal emission.
    """

    sl_path = os.environ['SL_CACHE_DIR']
    file_name = "%s/%s" % (sl_path, THERMAL_FNAME)
    f = open(file_name, 'r')
    lines_all = f.readlines()
    f.close()
    lines = [l for l in lines_all if l[0] != '#'] # allow any number of comment lines
    sep = ',' if ',' in lines[0] else ' '         # allow space- or comma-separated files

    thermal_wave = []
    thermal_list = []

    for aline in lines:
        (awave, athermal) = aline.strip().split(sep)
        awave = float(awave)
        athermal = float(athermal)
        thermal_wave.append(awave)
        thermal_list.append(athermal)

    # Set up the numpy arrays
    wave = np.array(thermal_wave, dtype='double')
    thermal_bg = np.array(thermal_list, dtype='double')
    return (wave, thermal_bg)
