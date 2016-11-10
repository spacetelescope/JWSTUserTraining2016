from pandeia.engine.helpers.sl_accessor_funcs import get_dateless_bg
from pandeia.engine.helpers.straylight_accessor import get_stray_light_bg
from pandeia.engine.helpers.thermal_accessor import get_thermal_bg
from pandeia.engine.helpers.bmg_accessor import get_in_field_bg

from pandeia.engine.custom_exceptions import BMGError
from pandeia.engine.custom_exceptions import BackgroundError
from pandeia.engine.custom_exceptions import DatelessBGError
from pandeia.engine.custom_exceptions import StraylightPositionError
from pandeia.engine.custom_exceptions import StraylightDataError

import pysynphot as pysyn
from operator import add
import numpy as np


'''
We have to resample the three background components because the infield waveset
is not the same as straylight and thermal.  So we normalize the waveset and then
resample each component.
'''


def bg_resample(merged, wave, flux):
    merged = np.array(merged)
    wave = np.array(wave)
    flux = np.array(flux)
    spec = pysyn.spectrum.ArraySourceSpectrum(wave=wave, flux=flux)
    f = np.ones(len(wave))
    filt = pysyn.spectrum.ArraySpectralElement(wave, f, waveunits='microns')
    obs = pysyn.observation.Observation(spec, filt, binset=merged, force='taper')
    flux = obs.binflux
    return list(flux)


def call_butler(ra, dec, date, level, ra_dec_str, date_str):

    if not date is None:
        # dated

        # stray light
        try:
            sl_wave, sl_bg = get_stray_light_bg(ra, dec, date, ra_dec_str, date_str)
        except (StraylightPositionError, StraylightDataError):
            raise
        except Exception as e:
            raise BackgroundError('Error calculating stray light background.', e)

        # thermal
        try:
            thermal_wave, thermal_bg = get_thermal_bg()
        except Exception as e:
            raise BackgroundError('Error calculating thermal background.', e)

        # in field
        try:
            if_wave, if_bg = get_in_field_bg(ra, dec, date, ra_dec_str, date_str)
        except BMGError as e:
            raise BackgroundError('BMGError calculating infield background: %s' % e)
        except Exception as e:
            raise BackgroundError('Error calculating infield background: %s' % e)

        # merge wavelengths
        # wave = astro_spectrum.merge_wavelengths(thermal_wave, if_wave)
        wave = if_wave

    else:

        # get all components via special get_dateless_bg route if no date given
        try:
            unused_iday, wave, if_bg, sl_bg, thermal_wave, thermal_bg = get_dateless_bg(ra, dec, level)
            # the get_dateless_bg() API returns iday (int day of yr) which we don't currently use
        except (DatelessBGError, StraylightPositionError, StraylightDataError):
            raise
        except Exception as e:
            raise BackgroundError('Error calculating dateless background.', e)
        # wave stands for both if_wave (infield) and sl_wave (straylight)

    # resample the flux for each onto the new set
    if not date is None:   # To Vicki: "if date is not None" also works but is ambiguous; might seem to mean "date is (not None)"
        sl_bg = bg_resample(wave, sl_wave, sl_bg)
    thermal_bg = bg_resample(wave, thermal_wave, thermal_bg)
    # if_bg = bg_resample(wave, if_wave, if_bg)

    # then sum the backgrounds
    combined_bg = map(add, thermal_bg, if_bg)
    combined_bg = map(add, combined_bg, sl_bg)

    # this will be written to a .npz file
    data_to_save = dict(
        straylight=sl_bg,
        thermal=thermal_bg,
        infield=if_bg,
        background=combined_bg,
        wavelength=wave
    )

    return [list(wave), combined_bg], data_to_save


def get_background(background):
    data_to_save = None

    # dated
    if background['bg_type'] == 'dated':
        ra = background['ra']
        dec = background['dec']
        ra_dec_str = background['ra_dec_str']
        date = background['date']
        date_str = background['date_str']
        background, data_to_save = call_butler(ra, dec, date, None, ra_dec_str, date_str)

    # dateless
    elif background['bg_type'] != None and background['bg_type'].lower() != 'none':
        ra = background['ra']
        dec = background['dec']
        ra_dec_str = background['ra_dec_str']
        level = background['bg_type'].upper()[0] # dateless code expects level of: L,M,H
        background, data_to_save = call_butler(ra, dec, None, level, ra_dec_str, "dateless")

    else:
        # Handle "none", but note that this is how we'd send positionless/dateless
        background = background['bg_type']

    return background, data_to_save
