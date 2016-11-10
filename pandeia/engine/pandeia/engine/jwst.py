# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import os

from astropy.io import fits

import numpy.ma as ma
import numpy as np

from .psf_library import PSFLibrary
from .telescope import Telescope
from .instrument import Instrument
from .io_utils import ref_data_interp, read_json
from .custom_exceptions import EngineInputError, RangeError, DataError


class JWST(Telescope):

    """
    This is currently a dummy class that is used for configuration file discovery. Could eventually
    contain JWST-specific methods.
    """
    pass


class JWSTInstrument(Instrument):

    """
    Generic JWST Instrument class
    """

    def __init__(self, mode=None, config={}, webapp=False, **kwargs):
        telescope = JWST()
        # these are the required sections and need to be passed via API in webapp mode
        self.instrument_pars = {}
        self.instrument_pars['detector'] = ["nexp", "ngroup", "nint", "readmode", "subarray"]
        self.instrument_pars['instrument'] = ["aperture", "disperser", "filter", "instrument", "mode"]
        self.api_parameters = list(self.instrument_pars.keys())

        # these are required for calculation, but ok to live with config file defaults
        self.api_ignore = ['dynamic_scene', 'max_scene_size', 'scene_size']

        Instrument.__init__(self, telescope=telescope, mode=mode, config=config, webapp=webapp, **kwargs)


class NIRSpec(JWSTInstrument):

    """
    Need to over-ride get_wave_range() for NIRSpec because the effective wavelength range
    depends on the blocking filter, the aperture, and the disperser.  Also need to overload __init__
    to handle some special MSA configuration needs.
    """

    def __init__(self, mode=None, config={}, webapp=False, **kwargs):
        JWSTInstrument.__init__(self, mode=mode, config=config, webapp=webapp, **kwargs)

        if self.mode == "msa":
            aperture = self.instrument['aperture']
            shutter_location = self.instrument['shutter_location']
            gap_config_file = self.aperture_config[aperture].pop('gap')
            gap_config = read_json(os.path.join(self.ref_dir, gap_config_file), raise_except=True)
            self.shutter_locations = gap_config.keys()
            try:
                self.aperture_config[aperture]["gap"] = gap_config[shutter_location]
            except KeyError as e:
                msg = "Shutter location not specified for MSA calculation: %s" % repr(e)
                raise DataError(value=msg)

    def get_wave_range(self):
        """
        Get the wavelength range of the current instrument configuration

        Returns
        -------
        range_dict: dict
            Contains the instrument wavelength range in microns described by:
                wmin - minimum wavelength
                wmax - maximum wavelength
        """
        disperser = self.instrument['disperser']
        aperture = self.instrument['aperture']
        filt = self.instrument['filter']
        gap = self.aperture_config[aperture]['gap']

        if disperser is not None:
            # get the wavelength range from the gap configuration
            if filt in gap[disperser]:
                # g140m and g140h have different gap configs for each blocking filter
                g_wmin = gap[disperser][filt]["wave_min"]
                g_wmax = gap[disperser][filt]["wave_max"]
            else:
                g_wmin = gap[disperser]["wave_min"]
                g_wmax = gap[disperser]["wave_max"]

            # get the wavelength range from the configuration file
            c_wmin = self.range[aperture][filt]["wmin"]
            c_wmax = self.range[aperture][filt]["wmax"]

            # get the wavelength range over which the disperser efficiency is known
            wave_blaze = self.get_wave_blaze()
            d_wmin = wave_blaze.min()
            d_wmax = wave_blaze.max()

            # get the wavelength range over which the filter throughput is known
            wave_filter = self.get_wave_filter()
            f_wmin = wave_filter.min()
            f_wmax = wave_filter.max()

            # compare filter and disperser wavelength ranges
            if f_wmax < d_wmin or d_wmax < f_wmin:
                raise RangeError(value="Disperser and Filter wavelength ranges do not overlap.")
            wmin = max(f_wmin, d_wmin, c_wmin, g_wmin)
            wmax = min(f_wmax, d_wmax, c_wmax, g_wmax)
        else:
            wmin = self.range[aperture]["wmin"]
            wmax = self.range[aperture]["wmax"]

        range_dict = {'wmin': wmin, 'wmax': wmax}
        return range_dict

    def create_gap_mask(self, wave):
        """
        Use the gap configuration and a wavelength vector, wave, to build a masked array that
        masks out the location of the gap.  Wavelengths between and including both gap endpoints
        will be masked out.

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to construct mask from

        Returns
        -------
        mask: numpy.ma 1D masked array
            1D array masked at the locations within the configured detector gap and 1.0 elsewhere
        """
        disperser = self.instrument['disperser']
        aperture = self.instrument['aperture']
        filt = self.instrument['filter']
        gap = self.aperture_config[aperture]['gap'][disperser]
        if filt in gap:
            gap_start = gap[filt]['gap_start']
            gap_end = gap[filt]['gap_end']
        else:
            gap_start = gap['gap_start']
            gap_end = gap['gap_end']

        if gap_start is not None and gap_end is not None:
            masked_wave = ma.masked_inside(wave, gap_start, gap_end)
            mask = masked_wave / wave
            mask = ma.filled(mask, 0.0)
        else:
            mask = 1.0
        return mask

    def get_internal_eff(self, wave):
        """
        Read in internal efficiency of NIRSpec. This is
        overloaded because the internal optical throughput is
        different for the NIRSpec IFU compared to MOS and Fixed Slit.

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate throughput onto

        Returns
        -------
        eff: numpy.ndarray or float
            Internal throughput as a function of wave
        """
        if self.mode == "ifu":
            eff = self._get_throughput(wave, 'internal_ifu')
        elif self.mode in ["msa", "fixed_slit"]:
            eff = self._get_throughput(wave, 'internal_msa')
        else:
            msg = "Internal efficiency not configured for NIRSpec mode %s." % self.mode
            raise EngineInputError(value=msg)

        return eff


class NIRCam(JWSTInstrument):

    def _loadpsfs(self):
        """
        For the bar-shaped coronagraphy masks we need to load the psf_library on a per-filter basis.
        The PSF files have the aperture and filter smooshed into one string so use that as the key.
        """
        if self.instrument['aperture'] in ('masklwb', 'maskswb'):
            psfs_key = "%s%s" % (self.instrument['aperture'], self.instrument['filter'])
        else:
            psfs_key = self.instrument['aperture']

        psf_path = os.path.join(self.ref_dir, "psfs")
        self.psf_library = PSFLibrary(path=psf_path, aperture=psfs_key)

    def get_filter_eff(self, wave):
        """
        over-ride the filter efficiency because the narrow-band filters are in the pupil wheel,
        and therefore also go through a broad-band filter in the filter wheel (which doesn't have a clear).

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate throughput onto

        Returns
        -------
        eff: numpy.ndarray or float
            Filter throughput as a function of wave
        """

        if not hasattr(self, 'double_filters'):
            msg = "NIRCam requires a mapping that describes which filters are actually a combination of two filters."
            raise DataError(value=msg)

        if self.instrument['filter'] in self.double_filters.keys():
            eff = self._get_throughput(wave, self.instrument['filter'])
            eff_pupil = self._get_throughput(wave, self.double_filters[self.instrument['filter']])
            eff *= eff_pupil
        else:
            eff = self._get_throughput(wave, self.instrument['filter'])
        return eff

    def get_detector_pars(self):
        """
        Need to over-ride get_detector_pars() to handle the two different detectors. Which one to use is keyed off
        of the configured aperture.
        """
        aperture = self.instrument['aperture']
        if aperture in ['lw', 'mask335r', 'mask430r', 'masklwb']:
            det_dict = self.det_pars['lw']
        elif aperture in ['sw', 'mask210r', 'maskswb']:
            det_dict = self.det_pars['sw']
        else:
            msg = "Unknown NIRCam aperture %s" % aperture
            raise EngineInputError(value=msg)

        return det_dict

    def get_internal_eff(self, wave):
        """
        Read in internal efficiency. For NIRCam there are separate internal efficiencies for the optics common
        to all modes, the throughput of the coronagrapher substrate, and the throughputs of the optical wedges
        that bring the coronagraphy elements into the field of view of the detectors.

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate throughput onto

        Returns
        -------
        eff: numpy.ndarray or float
            Internal throughput as a function of wave
        """
        base_eff = self._get_throughput(wave, 'internal')
        coronagraphy_eff = 1.0
        wedge_eff = 1.0
        dichroic_eff = 1.0

        # load the dichroic throughput based on which channel we're using
        ap = self.aperture_config[self.instrument['aperture']]
        if 'channel' not in ap:
            msg = "Channels not configured for NIRCam aperture %s." % self.instrument['aperture']
            raise DataError(value=msg)

        if ap['channel'] == 'sw':
            dichroic_eff = self._get_throughput(wave, 'dbs_sw')
        elif ap['channel'] == 'lw':
            dichroic_eff = self._get_throughput(wave, 'dbs_lw')
        else:
            msg = "Channel %s not configured for NIRCam aperture %s." % (ap['channel'], self.instrument['aperture'])
            raise DataError(value=msg)

        if self.instrument['mode'] == 'coronagraphy':
            coronagraphy_eff = self._get_throughput(wave, 'coronagraphy_substrate')
            if self.instrument['aperture'] in ('maskswb', 'mask210r'):
                wedge_eff = self._get_throughput(wave, 'sw_wedge_eff')
            if self.instrument['aperture'] == ('masklwb', 'mask335r', 'mask430r'):
                wedge_eff = self._get_throughput(wave, 'lw_wedge_eff')
        eff = base_eff * coronagraphy_eff * wedge_eff * dichroic_eff

        return eff

    def dispersion_axis(self):
        """
        The dispersion axis is either along rows (the X axis) or along columns (Y axis).  By default
        it is along the X axis.  However, the GRISMC grating disperses along the Y axis as a way to
        help mitigate source crowding and confusion.

        Returns
        -------
        disp_axis: str
            Axis along with spectra are dispersed.  Allowed values are 'x' or 'y'.
        """
        if self.projection_type == 'slitless' and self.instrument['disperser'] == "grismc":
            disp_axis = "y"
        else:
            disp_axis = "x"
        return disp_axis

    def bar_width(self, x):
        """
        Width of MASKLWB or MASKSWB in arcsec as a function of X. The width at the center of the FOV is taken from the
        configuration as a function of what filter is being used.

        Parameters
        ----------
        x: float
            X position (arcsec) in field of view

        Returns
        ------
        width: float
            Width of bar at X in arcsec
        """
        filt = self.instrument['filter']
        if filt not in self.bar_offsets:
            msg = "Invalid filter, %s, for MASKLWB/MASKSWB." % filt
            raise DataError(value=msg)
        # change the sign because webbpsf has X increase to left while we have it increasing to the right.
        # the bars narrow from left to right in the FOV.
        center = -1.0 * self.bar_offsets[filt]
        if self.instrument['aperture'] == 'maskswb':
            width = 0.2666 - 0.01777 * (center + x)
        elif self.instrument['aperture'] == 'masklwb':
            width = 0.5839 - 0.03893 * (center + x)
        else:
            msg = "bar_width() method only appropriate for MASKLWB and MASKSWB apertures."
            raise EngineInputError(value=msg)
        return width

    def get_detector_qe(self, wave):
        """
        Need to over-ride get_detector_qe() to handle the two different detectors. Which one to use is keyed off
        of the configured aperture.
        """
        aperture = self.instrument['aperture']
        try:
            channel = self.aperture_config[aperture]['channel']
        except KeyError as e:
            msg = "NIRCam aperture configuration must include which channel the aperture belongs to, sw or lw. (%s)" % e
            raise DataError(value=msg)
        if channel == 'lw':
            qe = self._get_throughput(wave, 'qe_lw')
        elif channel == 'sw':
            qe = self._get_throughput(wave, 'qe_sw')
        else:
            msg = "Unknown NIRCam channel %s" % channel
            raise DataError(value=msg)

        return qe


class NIRISS(JWSTInstrument):

    """
    Need to over-ride get_internal_eff() because there are two different wheels in NIRISS.
    The 'filter wheel' has a clear position with 100% transmission. The 'pupil wheel'
    (which somewhat confusingly also contains filters) has a clear position with some
    obstruction (the Pupil Alignment Reference - PAR). So if the active filter is in the
    'filter wheel', the transmission takes a hit from the PAR in the pupil wheel.

    Also need to override dispersion_axis() because the BG150C grating disperses along the
    Y axis vs. the X axis like every other JWST disperser.
    """

    def dispersion_axis(self):
        """
        The dispersion axis is either along rows (the X axis) or along columns (Y axis).  By default
        it is along the X axis.  However, the GR150C grating disperses along the Y axis as a way to
        help mitigate source crowding and confusion.

        Returns
        -------
        disp_axis: str
            Axis along with spectra are dispersed.  Allowed values are 'x' or 'y'.
        """
        if self.instrument['disperser'] in ("gr150c", "gr700xd"):
            disp_axis = "y"
        else:
            disp_axis = "x"
        return disp_axis

    def get_internal_eff(self, wave):
        """
        Calculate NIRISS internal efficiency

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate efficiency onto

        Returns
        -------
        eff: numpy.ndarray or float
            Internal efficiency as a function of wave
        """
        filt = self.instrument['filter']
        if self.mode in ['ami', 'imaging']:
            # Pupil wheel
            if filt in ["f090w", "f115w", "f158m", "f140m", "f150w", "f200w"]:
                pupil_file = os.path.join(self.ref_dir, self.paths["clear"])
            # Filter wheel
            elif filt in ["f277w", "f444w", "f356w", "f430m", "f380m", "f480m"]:
                pupil_file = os.path.join(self.ref_dir, self.paths["clearp"])

            pupil_eff = ref_data_interp(pupil_file, wave, 'throughput')
        else:
            pupil_eff = 1.0

        # handle NRM transmission
        if self.mode == 'ami':
            nrm_file = os.path.join(self.ref_dir, self.paths["nrm"])
            nrm_eff = ref_data_interp(nrm_file, wave, 'throughput')
        else:
            nrm_eff = 1.0

        internal_file = os.path.join(self.ref_dir, self.paths["internal"])
        internal_eff = ref_data_interp(internal_file, wave, 'throughput')

        # Total internal efficiency are the internal reflections and the pupil/filter wheel transmissions
        eff = internal_eff * pupil_eff * nrm_eff

        return eff

    def get_extraction_mask(self, order):
        """
        Each SOSS order has its own extraction mask. Use the specified order to build the key to look up the
        FITS file containing the mask.

        Parameters
        ----------
        order: int
            Order whose mask to read

        Returns
        -------
        mask: 2D np.ndarray
            Mask data
        """
        if order not in (1, 2, 3):
            msg = "SOSS order %d is not valid." % order
            raise EngineInputError(value=msg)

        key = "gr700xd_%d_mask" % order

        # substrip96 is a special case where there's only one possible mask
        if self.detector['subarray'] == 'substrip96':
            key += "96"

        path = self.paths.get(key, None)

        if path is None:
            msg = "No mask configured for SOSS order %d." % order
            raise DataError(value=msg)

        mask_file = os.path.join(self.ref_dir, path)
        try:
            mask = fits.getdata(mask_file)
        except Exception as e:
            msg = "Error reading mask data for SOSS order %d: %s" % (order, repr(e))
            raise DataError(value=msg)

        return mask

    def get_trace(self, wave):
        """
        Read in spectral trace offsets from center of FOV. Currently wavelength-dependent spatial distortion is
        only required for SOSS mode. Other modes simply return 0's.

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate/trim trace data onto

        Returns
        -------
        trace: numpy.ndarray or float
            Spectral trace offset from center of FOV as a function of wave
        """
        if self.mode == "soss":
            disperser = self._get_disperser()
            key = "%s_wavepix" % disperser
            # handle the special case of substrip96
            if self.detector['subarray'] == 'substrip96':
                key += "96"

            # SOSS modes requires trace reference data to work properly. so raise exception if we can't load it.
            try:
                trace = self._interp_refdata(wave, key, colname='trace')
            except DataError as e:
                msg = "Spectral trace data missing for NIRISS SOSS, %s."
                raise DataError(value=msg)
        else:
            trace = np.zeros_like(wave)
        return trace

    def get_detector_pixels(self, wave):
        """
        Read in detector pixel positions for each wavelength in wave_pix

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate/trim pixel position data onto

        Returns
        -------
        pixels: numpy.ndarray or None
            Detector pixel positions corresponding to each element of wave
        """
        disperser = self._get_disperser()
        if disperser is not None:
            key = "%s_wavepix" % disperser
            # handle the special case of substrip96
            if self.detector['subarray'] == 'substrip96':
                key += "96"
            try:
                pixels = self._interp_refdata(wave, key, colname='detector_pixels', default=np.arange(len(wave)))
            except DataError as e:
                pixels = None
        else:
            pixels = None
        return pixels


class MIRI(JWSTInstrument):

    """
    The MIRI internal efficiency, detector readmodes, etc. are more complex, and different than the other instruments,
    so some methods are redefined.
    """

    def get_detector_pars(self):
        """
        Get detector parameters which are keyed by short wavelength (sw), long wavelength (lw),
        and imager.  Use configured disperser and aperture to figure out which to use.

        Returns
        -------
        det_dict: dict
            rn - readout noise (float)
            ipc - detector has inter-pixel capacitance (bool)
            dark_current - dark current (float)
            fullwell - pixel capacity in electrons (float)
            rn_correlation - detector has correlated readout noise (bool)
            ff_electrons - number of counts in flat-field frame (float)
        """
        disperser = self.instrument["disperser"]
        aperture = self.instrument["aperture"]

        if aperture == 'ch1' or aperture == 'ch2':
            det_dict = self.det_pars["sw"]
        elif aperture == 'ch3' or aperture == 'ch4':
            det_dict = self.det_pars["lw"]
        elif disperser is None or disperser == 'p750l':
            det_dict = self.det_pars["imager"]

        return det_dict

    def get_detector_qe(self, wave):
        """
        Read in detector quantum efficiency. Need to overload this here to know when to use
        "sw", "lw", or "imager".

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate throughput onto

        Returns
        -------
        qe: numpy.ndarray or float
            Quantum efficiency as a function of wave
        """
        aperture = self.instrument["aperture"]
        if aperture in ['ch1', 'ch2']:
            key = "miri_sw_qe"
        elif aperture in ['ch3', 'ch4']:
            key = "miri_lw_qe"
        else:
            key = "miri_imager_qe"

        qe = self._get_throughput(wave, key)
        return qe

    def get_dichroic(self, type, level):
        """
        MIRI dichroic reference files are keyed by type (refl or trans), level, and disperser

        Parameters
        ----------
        type: str
            Either 'refl' for reflective or 'trans' for transmitting
        level: int
            Level 1, 2, or 3

        Returns
        -------
        key: str
            Key pointing to the dichroic efficiency data in self.paths
        """
        disperser = self.instrument['disperser']
        if disperser == 'short':
            mrs_setting = 'a'
        elif disperser == 'medium':
            mrs_setting = 'b'
        elif disperser == 'long':
            mrs_setting = 'c'

        key = "sd" + str(level) + mrs_setting + "_" + type
        return key

    def get_internal_eff(self, wave):
        """
        Calculate MIRI internal efficiency which is rather more complicated than the other instruments

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate efficiency onto

        Returns
        -------
        eff: numpy.ndarray or float
            Internal efficiency as a function of wave
        """
        aperture = self.instrument['aperture']
        mirror_eff = self.mirror_eff
        mirror_cont = self.mirror_cont
        n_refl = self.n_reflections[aperture]
        if self.mode == 'mrs':
            if aperture == 'ch1':
                d1r_key = self.get_dichroic('refl', 1)
                d1r_eff = self._get_throughput(wave, d1r_key)
                refl_eff = mirror_eff ** n_refl
                internal_eff = d1r_eff * refl_eff

            elif aperture == 'ch2':
                d1t_key = self.get_dichroic('trans', 1)
                d2r_key = self.get_dichroic('refl', 2)
                d1t_eff = self._get_throughput(wave, d1t_key)
                d2r_eff = self._get_throughput(wave, d2r_key)
                refl_eff = mirror_eff ** n_refl
                internal_eff = d1t_eff * d2r_eff * refl_eff

            elif aperture == 'ch3':
                d1t_key = self.get_dichroic('trans', 1)
                d2t_key = self.get_dichroic('trans', 2)
                d3r_key = self.get_dichroic('refl', 3)
                d1t_eff = self._get_throughput(wave, d1t_key)
                d2t_eff = self._get_throughput(wave, d2t_key)
                d3r_eff = self._get_throughput(wave, d3r_key)
                refl_eff = mirror_eff ** n_refl
                internal_eff = d1t_eff * d2t_eff * d3r_eff * refl_eff

            elif aperture == 'ch4':
                d1t_key = self.get_dichroic('trans', 1)
                d2t_key = self.get_dichroic('trans', 2)
                d3t_key = self.get_dichroic('trans', 3)
                d1t_eff = self._get_throughput(wave, d1t_key)
                d2t_eff = self._get_throughput(wave, d2t_key)
                d3t_eff = self._get_throughput(wave, d3t_key)
                refl_eff = mirror_eff ** n_refl
                internal_eff = d1t_eff * d2t_eff * d3t_eff * refl_eff

            else:
                raise EngineInputError(value='Invalid aperture for MIRI: %s' % aperture)

        else:
            refl_eff = mirror_eff ** n_refl
            internal_eff = wave * 0.0 + refl_eff

        # mirror contamination factor
        internal_eff = internal_eff * mirror_cont

        return internal_eff

    def get_disperser_eff(self, wave):
        """
        Overloaded here because disperser efficiency is keyed off of the aperture rather than disperser

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate throughput onto

        Returns
        -------
        disperser_eff: numpy.ndarray or float
            Disperser efficiency as a function of wave
        """
        if self.instrument['disperser'] is not None:
            key = self.instrument['aperture']
            disperser_eff = self._get_throughput(wave, key)
        else:
            disperser_eff = 1.
        return disperser_eff

    def _get_dispersion_key(self):
        """
        Overload this because the key is constructed from both the aperture and disperser rather
        than the disperser alone.

        Returns
        -------
        key: str
            Key used to get dispersion file out of self.paths
        """
        disperser = self.instrument['disperser']
        aperture = self.instrument['aperture']
        key = "%s_%s_disp" % (aperture, disperser)
        return key
