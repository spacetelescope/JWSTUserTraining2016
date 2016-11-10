# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import os
import numpy as np
from astropy.io import ascii, fits
from astropy.table import unique

from .config import DefaultConfig
from .custom_exceptions import DataConfigurationError, EngineInputError
from .utils import get_key_list, get_dict_from_keys, recursive_subclasses, merge_data
from .io_utils import read_json
from . import config as cf

default_refdata_directory = cf.default_refdata_directory


class Extinction(DefaultConfig):

    """
    Class for handling configuration and application of spectrum NormalizationFactory

    Attributes
    ----------
    law: string
        Specifies which extinction law to use
    unit: string
        Units of extinction, either 'nh' for hydrogen column density (in cm^-2) or 'mag' for magnitudes of
        extinction in a provided bandpass, 'bandpass'
    value: float
        Amount of extinction in unit
    bandpass: string
        If unit is 'mag', this is provided to specify the bandpass over which the extinction is integrated
    extinction_wave: np.array
        Wavelength vector of the configured extinction law (microns)
    extinction_curve: np.array
        Configured extinction law in units of cross-section per hydrogen atom, cm^2/H
    """

    def __init__(self, webapp=False, config={}, **kwargs):
        """
        How an extinction law is configured depends on the defaults, the per-law defaults, and any input parameters.

        Parameters
        ----------
        webapp: bool
            Determines whether strict API checks should be performed
        config: dict
            Extra configuration information in engine input API dict format
        **kwargs: keyword/value pairs
            Extra configuration information in kwargs format
        """
        self.webapp = webapp

        # send webapp=False here since we do the API checking further below
        DefaultConfig.__init__(self, webapp=False, config=config, **kwargs)

        all_config = merge_data(config, dict(**kwargs))

        if not hasattr(self, "law"):
            self.law = self.__class__.__name__.lower()

        if self.law in self.laws:
            law_conf = self.laws[self.law]
            if "defaults" in law_conf:
                law_defaults = law_conf.pop("defaults")
            else:
                law_defaults = {}
            law_conf = merge_data(law_conf, law_defaults, config, dict(**kwargs))
            self.__dict__ = merge_data(self.__dict__, law_conf)
        else:
            msg = "Unsupported extinction law, %s" % self.law
            raise EngineInputError(value=msg)

        # do some sanity checks to make sure inputs make sense
        try:
            self._sanity_checks()
        except AttributeError as e:
            self.warnings['no_sanity_check'] = warning_messages['no_sanity_check'] % (self.__class__.__name__, e)

        # add 'bandpass' to api_parameters if we need it
        if self.unit == "mag":
            self.api_parameters.append('bandpass')

        # we need to run the API checks here after merging in the normalization type-specific defaults
        if webapp:
            self._api_checks(all_config)

        # load the extinction curve, i.e. extinction cross section versus wavelength
        self.extinction_wave, self.extinction_curve = self._load_curve()

    def _sanity_checks(self):
        # make sure we're using a valid unit
        if self.unit not in self.units:
            msg = "Invalid extinction unit, %s, specified ." % self.unit
            raise EngineInputError(value=msg)

        # make sure if we're using a bandpass, we have it defined
        if self.unit == "mag":
            if self.bandpass not in self.bandpasses:
                msg = "Invalid extinction bandpass, %s, specified ." % self.bandpass
                raise EngineInputError(value=msg)

    def _get_config(self):
        """
        Read default configuration from JSON

        Returns
        -------
        config: dict
            All desired class attributes should have defaults defined in the config file
        """
        # get the configuration data from the pandeia reference data directory
        ref_dir = os.path.join(default_refdata_directory, "extinction")
        config = read_json(os.path.join(ref_dir, "config.json"), raise_except=True)

        # pop the defaults entry out
        if "defaults" in config:
            defaults = config.pop("defaults")
            config.update(defaults)
        else:
            msg = "No extinction defaults defined."
            raise DataConfigurationError(value=msg)

        return config

    def _load_curve_data(self):
        """
        Use the configuration data for a given curve to configure the astropy.io.ascii fixed-width reader and load the data.

        Returns
        -------
        data: dict-like
            Data columns read from extinction data file.
        """
        curve_config = self.laws[self.law]
        filename = os.path.join(default_refdata_directory, "extinction", "curves", curve_config['filename'])
        data = ascii.read(
            filename,
            Reader=ascii.FixedWidth,
            data_start=curve_config['data_start'],
            names=curve_config['names'],
            col_starts=curve_config['col_starts'],
            col_ends=curve_config['col_ends'],
            guess=False
        )
        # make sure the data is sorted in order of increasing wavelength. 'data' is an astropy.table.Table object so has
        # built-in method to do this.
        data.sort(['wave'])

        # some of the files have duplicate rows so use astropy.table.unique() to weed them out, too
        uniq_data = unique(data, keys='wave')
        return uniq_data

    def extinction(self, wave, flux):
        """
        This is to apply the actual extinction to an input (wave, flux) set of vectors.  This is primarily
        based on comments in https://github.com/STScI-SSB/pandeia/issues/503#issuecomment-232136894.

        Parameters
        ----------
        wave - vector of wavelengths in ?? units
        flux - vector of flux (same shape as wave vector)

        Returns
        -------
        wave - a copy of the input vector of wavelengths
        flux - extincted version of the flux.

        """
        # interpolate the extinction curve onto the input wavelength set
        extinction_curve_interp = np.interp(wave, self.extinction_wave, self.extinction_curve)

        if self.unit == "nh":
            # if self.unit is 'nh', then self.value is a hydrogen column density in cm^-2 and we use it to scale the extinction
            # curve, C_ext, directly.
            nh_cext = self.value * extinction_curve_interp
        elif self.unit == "mag":
            # the other possible unit is 'mag' in which case we use the provided bandpass to
            # normalize the extinction curve. self._sanity_checks will raise exception before we reach
            # this point if self.unit is not 'nh' or 'mag'.
            keys = get_key_list(self.bandpass, separator=',')  # should be [photsys, filter]
            bp_filename = get_dict_from_keys(self.bandpasses, keys)['filename']

            bp_path = os.path.join(default_refdata_directory, "normalization", "bandpass", bp_filename)

            bp_fits = fits.open(bp_path)
            bp_wave, bp_throughput = [], []
            for r in bp_fits[1].data:
                bp_wave.append(r[0])
                bp_throughput.append(r[1])
            bp_wave, bp_throughput = np.array(bp_wave), np.array(bp_throughput)
            bp_wave = bp_wave / 10000.0  # convert to microns

            # Need to interpolate the bandpass / throughput bp.wave and bp.throughput onto the self.wave
            # values and then we can do the normalization.
            # np.interp is piecewise linear
            bp_throughput_interp = np.interp(self.extinction_wave, bp_wave, bp_throughput)

            # Calculate the bandpass normalized extinction curve value
            #
            # There's more work to be done here. The normalization formally also requires the SED of the source as well.
            # This simplification is good enough for now and matches what's most commonly used.  However, for complicated
            # SEDs subject to large amounts of complicated extinction, the more formal approach will be required.
            # The work required for this is described in https://github.com/STScI-SSB/pandeia/issues/1884.
            bandpass_normalized_c_ext = sum(np.multiply(bp_throughput_interp, self.extinction_curve)) / \
                                        sum(bp_throughput_interp)

            # Calculate A_lambda
            # see https://github.com/STScI-SSB/pandeia/issues/503#issuecomment-232136894
            A_lambda = self.value * extinction_curve_interp / bandpass_normalized_c_ext
            # convert from A_lambda in magnitudes to NH * C_ext
            nh_cext = A_lambda / 1.08574
        else:
            msg = "Unsupported extinction unit, %s" % self.unit
            raise EngineInputError(value=msg)

        # Apply the extinction to the input flux
        extincted_flux = np.exp(-nh_cext) * flux

        return wave, extincted_flux


class WD2001(Extinction):

    """
    This class implements the extinction models of Weingartner & Draine 2001, ApJ, 548, 296 with updates where applicable
    from Draine 2003a, b, c.  The extinction laws from these models are tabulated and cover wavelengths from 1 A to 1 cm
    (1e-4 to 1e4 microns).  These models calculate the interstellar extinction using a dust model of carbonaceous grains
    and amorphous silicate grains.  The carbonaceous grains are like PAHs when small and like graphite when large.  The
    model has been calculated for different grain size distributions to match the observed extinction along a variety of
    lines of sight.  The extinction laws calculated from this model include:
        - Milky Way for an R_V value of 3.1 (default)
        - Milky Way for an R_V value of 4.0
        - Milky Way for an R_V value of 5.5
        - High-latitude molecular cloud HD 210121 with C/H = b_C = 40 ppm in log-normal size dists
        - Average for the LMC with C/H = b_C = 20 ppm in log-normal size dists
        - LMC with C/H = b_C = 10 ppm in log-normal size dists
        - SMC bar with C/H = b_C = 0 ppm in log-normal size dists
    """

    def _load_curve(self):
        """
        The WD2001 data files provide the extinction cross section directly in units of cm^2/H. No further
        conversion is needed except for checking of the wavelength order.

        Returns
        -------
        wave: np.array
            Wavelength vector of extinction curve
        c_ext: np.array
            Extinction curve
        """
        data = self._load_curve_data()

        wave = np.asarray(data['wave'])
        c_ext = np.asarray(data['c_ext'])

        return wave, c_ext


class Chapman09(Extinction):

    """
    This class implements the extinction law model of Chapman et al. 2009, ApJ, 690, 496.  It is most appropriate for dense
    molecular clouds and has been derived using mid-IR observations of a set of three molecular clouds: Ophiuchus, Perseus,
    and Serpens.
    """

    def _load_curve(self):
        """
        The Chapman et al. (2009) model provides the absorption and scattering cross sections separately. These need to be
        combined and then converted from cm^2/g to cm^2/H.

        Returns
        -------
        wave: np.array
            Wavelength vector of extinction curve
        c_ext: np.array
            Extinction curve
        """
        data = self._load_curve_data()

        # get the wave, C_absorption, and C_scattering data out
        wave = np.asarray(data['wave'])
        ca = np.asarray(data['ca'])
        cs = np.asarray(data['cs'])

        # add ca and cs to get total cross section and convert cm^2/g to cm^2/H.  use the WD01 R_V=5.5 factor of 2.2e-26 g/H
        # times the dust/ice mass fraction of 1.5 => 3.3e-26 g/H
        c_ext = 3.3e-26 * (ca + cs)

        return wave, c_ext


class MWRV31(WD2001):
    """
    This class wraps the md_rv_31 extinction law that is derived from the Weingartner & Draine model for Milky Way
    extinction with R_V = 3.1
    """
    pass


class MWRV40(WD2001):
    """
    This class wraps the md_rv_40 extinction law that is derived from the Weingartner & Draine model for Milky Way
    extinction with R_V = 4.0
    """
    pass


class MWRV55(WD2001):
    """
    This class wraps the md_rv_55 extinction law that is derived from the Weingartner & Draine model for Milky Way
    extinction with R_V = 5.5
    """
    pass


class HD210121(WD2001):
    """
    This class wraps the hd210121 extinction law that is derived from the Weingartner & Draine model for extinction along
    the line of sight to the high latitude molecular cloud HD 210121.
    """
    pass


class LMCavg(WD2001):
    """
    This class wraps the lmc_avg extinction law that is derived from the Weingartner & Draine model for average extinction
    within the LMC.
    """
    pass


class LMC2(WD2001):
    """
    This class wraps the lmc_2 extinction law that is derived from the Weingartner & Draine model for a modified versions
    of the LMC extinction model with C/H = b_C = 10 ppm in log-normal size distributions.
    """
    pass


class SMCbar(WD2001):
    """
    This class wraps the smc_bar extinction law that is derived from the Weingartner & Draine model for extinction within
    the bar region of the SMC.
    """
    pass


def ExtinctionFactory(webapp=False, config={}, **kwargs):

    """
    Function to take extinction configuration data and build/return a configured instance of the appropriate
    subclass of Extinction.

    Parameters
    ----------
    webapp: bool
        Switch to toggle strict API checking
    config: dict
        Configuration data in engine API dict format
    **kwargs: keyword/value pairs
        Additional configuration data
    """
    all_config = merge_data(config, dict(**kwargs))
    if 'law' in all_config:
        method = all_config['law'].replace('_', '')
    else:
        msg = "Must specify extinction law to use."
        raise EngineInputError(value=msg)

    laws = recursive_subclasses(Extinction)
    methods = [l.__name__.lower() for l in laws]
    law_map = dict(list(zip(methods, laws)))

    if method not in methods:
        msg = "Unsupported or not yet implemented extinction law: %s" % method
        raise EngineInputError(value=msg)
    else:
        cls = law_map[method](webapp=webapp, config=config, **kwargs)
        return cls
