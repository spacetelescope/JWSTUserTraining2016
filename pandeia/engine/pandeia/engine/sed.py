# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import os
import numpy as np

from astropy import units as u
from astropy.analytic_functions import blackbody_nu

import pysynphot as psyn

from .config import DefaultConfig
from .custom_exceptions import DataConfigurationError, EngineInputError
from .pandeia_warnings import sed_warning_messages as warning_messages
from .utils import recursive_subclasses, merge_data
from .constants import pandeia_waveunits, pandeia_fluxunits
from .io_utils import read_psyn_spectrum
from . import config as cf
from pandeia.engine.io_utils import read_json

default_refdata_directory = cf.default_refdata_directory


class SED(DefaultConfig):

    """
    Base class for configuration of spectral energy distributions
    """

    def __init__(self, webapp=False, config={}, **kwargs):
        """
        How an SED is configured depends on the defaults and any input configuration.

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

        # API checks are handled below so don't do it here...
        DefaultConfig.__init__(self, webapp=False, config=config, **kwargs)

        if not hasattr(self, "sed_type"):
            self.sed_type = self.__class__.__name__.lower()

        if self.sed_type in self.sed_types:
            sed_config = self.sed_types[self.sed_type]
            if "defaults" in sed_config:
                type_defaults = sed_config.pop("defaults")
                sed_config = merge_data(sed_config, type_defaults)
            else:
                msg = "No defaults defined for SED format %s" % self.sed_type
                raise DataConfigurationError(value=msg)

        self.__dict__ = merge_data(sed_config, self.__dict__)

        # do some sanity checks to make sure inputs make sense
        try:
            self._sanity_checks()
        except AttributeError as e:
            self.warnings['no_sanity_check'] = warning_messages['no_sanity_check'] % (self.__class__.__name__, e)

        # we need to run the API checks here after merging in the SED type-specific defaults
        if webapp:
            all_config = merge_data(config, dict(**kwargs))
            self._api_checks(all_config)

    def _get_config(self):
        """
        Read default configuration from JSON

        Returns
        -------
        config: dict
            All desired class attributes should have defaults defined in the config file
        """
        # use this trick to key the configuration file name off the name of the instantiated subclass
        ref_dir = os.path.join(default_refdata_directory, "sed")
        config = read_json(os.path.join(ref_dir, "config.json"), raise_except=True)
        defaults = config.pop("defaults")
        config.update(defaults)
        # nuke some unused things that are there for UI use
        for k in ["sed_families"]:
            del config[k]
        return config


class Input(SED):

    """
    Base class for SEDs provided as arrays
    """

    def __init__(self, webapp=False, config={}, **kwargs):
        """
        Parameters
        ----------
        webapp: bool
            Determines whether strict API checks should be performed
        config: dict
            Extra configuration information in engine input API dict format
        **kwargs: keyword/value pairs
            Extra configuration information in kwargs format
        """
        SED.__init__(self, webapp=webapp, config=config, **kwargs)
        self.wave, self.flux = self.spectrum
        if isinstance(self.wave, (list, tuple)):
            self.wave = np.array(self.wave)
        if isinstance(self.flux, (list, tuple)):
            self.flux = np.array(self.flux)

        if np.any(np.diff(self.wave) < 0):
            indices = np.where(np.diff(self.wave) < 0)[0]
            msg = "Input spectrum must be sorted by wavelength. Out-of-order indices: %s" % repr(indices)
            raise EngineInputError(value=msg)

        if np.any(np.diff(self.wave) == 0.0):
            indices = np.where(np.diff(self.wave) == 0.0)[0]
            msg = "Input spectrum contains duplicate wavelengths. Duplicate indices: %s" % repr(indices)
            raise EngineInputError(value=msg)

        self.wmin = self.wave.min()
        self.wmax = self.wave.max()

    def _sanity_checks(self):
        """
        Make sure we're given a spectrum and it has elements for wavelength and flux.
        We don't do any unit handling (yet) so units are assumed to be pandeia_waveunits
        and pandeia_fluxunits
        """
        if not hasattr(self, "spectrum"):
            msg = "Missing input spectrum."
            raise DataConfigurationError(value=msg)

        if len(self.spectrum) != 2:
            msg = "Must provide input spectrum as a 2 element list-like where 0th index is wavelength and 1st is flux."
            raise EngineInputError(value=msg)


class Analytic(SED):

    """
    Base class for analytic SEDs
    """

    def __init__(self, webapp=False, config={}, **kwargs):
        """
        Parameters
        ----------
        webapp: bool
            Determines whether strict API checks should be performed
        z: float
            Used to adjust wmin/wmax to account for source redshift.
        config: dict
            Extra configuration information in engine input API dict format
        **kwargs: keyword/value pairs
            Extra configuration information in kwargs format
        """
        SED.__init__(self, webapp=webapp, config=config, **kwargs)
        # the wmin and wmax configuration values are in the JWST frame. we need to back out
        # the redshift to calculate the analytic spectra in the correct rest frame.
        self.wmin /= 1.0 + self.z
        self.wmax /= 1.0 + self.z
        self.wave = self.create_wave()
        self.flux = self.create_flux()

    def _sanity_checks(self):
        """
        Check various configuration items to make sure they make sense and raise exception if not.
        """
        if hasattr(self, "unit") and hasattr(self, "units"):
            if self.unit not in self.units:
                msg = "Unsupported unit, %s, for analytic spectrum type %s" % (self.unit, self.__class__.__name__)
                raise EngineInputError(value=msg)

        for par in ["wmin", "wmax", "sampling", "z"]:
            if not hasattr(self, par):
                msg = "%s spectrum missing %s." % (self.__class__.__name__, par)
                raise DataConfigurationError(value=msg)

    def create_wave(self):
        # Create a wavelength grid with constant lambda/delta_lambda
        try:
            k = -self.sampling * np.log(self.wmin)
            nw = int(self.sampling * np.log(self.wmax) + k)
            wave = np.exp((np.arange(nw) - k) / self.sampling)
        except ValueError as e:
            msg = "Error creating wavelength grid for wmin=%f, wmax=%s with sampling=%d: %s" % (
                self.wmin,
                self.wmax,
                self.sampling,
                e
            )
            raise EngineInputError(value=msg)
        return wave

    def create_flux(self):
        flux = np.zeros(len(self.wave))
        return flux


class Powerlaw(Analytic):

    """
    Class implementing a power-law continuum of the form Flux ~ lam ** index
    """

    def _sanity_checks(self):
        """
        Run sanity checks specify to power-law SEDs
        """
        for par in ['unit', 'units', 'index']:
            if not hasattr(self, par):
                msg = "Missing required parameter %s for %s." % (par, self.__class__.__name__)
                raise DataConfigurationError(value=msg)

    def create_flux(self):
        # effectively normalize to 1 mJy at 1 micron and use Normalize() to scale from there
        flux = self.wave ** self.index
        if self.unit == 'flam':
            # pysynphot defaults to photlam which is scaled flam
            sp = psyn.ArraySpectrum(self.wave, flux)
            sp.convert(pandeia_waveunits)
            sp.convert(pandeia_fluxunits)
            flux = sp.flux

        return flux


class Flat(Analytic):

    """
    Class implementing a flat continuum in either Fnu or Flam
    """

    def _sanity_checks(self):
        """
        Run sanity checks specify to Flat SEDs
        """
        Analytic._sanity_checks(self)
        for par in ['unit', 'units']:
            if not hasattr(self, par):
                msg = "Missing required parameter %s for %s." % (par, self.__class__.__name__)
                raise DataConfigurationError(value=msg)

        if self.unit not in self.units:
            msg = "Unsupported unit, %s, for %s." % (par, self.__class__.__name__)
            raise EngineInputError(value=msg)

    def create_flux(self):
        flux = np.ones(self.wave.size)
        if self.unit == 'flam':
            # pysynphot defaults to photlam which is scaled flam
            sp = psyn.ArraySpectrum(self.wave, flux)
            sp.convert(pandeia_waveunits)
            sp.convert(pandeia_fluxunits)
            flux = sp.flux

        return flux


class BlackBody(Analytic):

    """
    Class implementing a Plank blackbody continuum
    """

    def _sanity_checks(self):
        """
        Run sanity checks specify to blackbody SEDs
        """
        Analytic._sanity_checks(self)
        if not hasattr(self, 'temp'):
            msg = "Missing required parameter 'temp' for %s." % self.__class__.__name__
            raise DataConfigurationError(value=msg)

    def create_flux(self):
        # use astropy analytic function and astropy.units to specify microns and K.
        # flux is in erg cm-2 s-1 hz-1 sr-1, but will be scaled by normalization.
        flux = blackbody_nu(self.wave * u.micron, self.temp * u.K)
        return flux.value


class No_Continuum(Analytic):

    """
    Class implementing zero continuum (for use with emission lines).
    Implement as just an alias for the base class.
    """

    pass


class Parameterized(SED):

    """
    Class implementing SED's that are calculated from a model using a set of input parameters.
    The parameters can either come from a specfied set contained within a configured catalog
    (usual for webapp mode) or as input configuration data.
    """

    def __init__(self, webapp=False, config={}, **kwargs):
        """
        How an SED is configured depends on the defaults and any input configuration.

        Parameters
        ----------
        webapp: bool
            Determines whether strict API checks should be performed
        config: dict
            Extra configuration information in engine input API dict format
        **kwargs: keyword/value pairs
            Extra configuration information in kwargs format
        """
        SED.__init__(self, webapp=webapp, config=config, **kwargs)
        spectra_file = os.path.join(default_refdata_directory, "sed", self.spectra["config"])
        self.spectra = read_json(spectra_file)
        self.wave, self.flux = self.get_spectrum()
        self.wmin = self.wave.min()
        self.wmax = self.wave.max()

    def _sanity_checks(self):
        """
        If self.webapp=True, make sure a key is provided and in the configured catalog.  If self.webapp=False,
        make sure the required parameters are
        """
        if self.webapp:
            if not hasattr(self, "key"):
                msg = "Must provide a key for looking up a SED."
                raise DataConfigurationError(value=msg)
        else:
            for p in self.api_parameters:
                if not hasattr(self, p):
                    msg = "Missing parameter, %s, for %s" % (p, self.__class__.__name__)
                    raise EngineInputError(value=msg)


class Phoenix(Parameterized):

    """
    Class implementing the Phoenix parameterized stellar spectrum models. The inputs to the models
    are:

    teff - effective temperature (K; range 2000 to 70000)
    metallicity - log10 of the metallicity in units of solar metallicity (range -4.0 to 0.5)
    log_g - surface gravity (log10 in cgs units; range 0.0 to 5.5)

    If used with webapp=True, the parameters will come from a configured look-up table
    using a provided key. Otherwise the model parameters will be taken from default configuration or input.
    """

    def __init__(self, webapp=False, config={}, **kwargs):
        if webapp:
            self.api_parameters = ["sed_type", "key"]
            self.api_ignore = ["spectrum_id", "teff", "metallicity", "log_g", "comment", "display_string", "spectra", "z"]
        else:
            self.api_parameters = ["sed_type", "teff", "metallicity", "log_g"]
            self.api_ignore = ["spectrum_id", "key", "comment", "display_string", "spectra", "z"]

        Parameterized.__init__(self, webapp=webapp, config=config, **kwargs)

    def get_spectrum(self):
        """
        Use pysynphot.Icat() to calculate spectrum and convert to microns/mJy.  If self.webapp=True,
        use key to look up a specified set of parameters.  Otherwise, get them from the configured
        attributes.
        """
        if self.webapp:
            if self.key not in self.spectra:
                msg = "Provided SED key, %s, not supported." % self.key
                raise EngineInputError(value=msg)
            pars = self.spectra[self.key]
            m = pars['metallicity']
            t = pars['teff']
            g = pars['log_g']
        else:
            m = self.metallicity
            t = self.teff
            g = self.log_g
        sp = psyn.Icat("phoenix", t, m, g)
        sp.convert("microns")
        sp.convert("mjy")
        return sp.wave, sp.flux


class Mapped(SED):

    """
    Class implementing SED's that are looked-up via an input 'key' and loaded from a file.
    """

    def __init__(self, webapp=False, config={}, **kwargs):
        """
        How an SED is configured depends on the defaults and any input configuration.

        Parameters
        ----------
        webapp: bool
            Determines whether strict API checks should be performed
        config: dict
            Extra configuration information in engine input API dict format
        **kwargs: keyword/value pairs
            Extra configuration information in kwargs format
        """
        SED.__init__(self, webapp=webapp, config=config, **kwargs)
        spectra_file = os.path.join(default_refdata_directory, "sed", self.spectra["config"])
        self.spectra = read_json(spectra_file)
        self.wave, self.flux = self.get_spectrum()
        self.wmin = self.wave.min()
        self.wmax = self.wave.max()

    def _sanity_checks(self):
        """
        Make sure a key is provided
        """
        if not hasattr(self, "key"):
            msg = "Must provide a key for looking up a SED."
            raise DataConfigurationError(value=msg)

    def get_spectrum(self):
        """
        Use self.key to grab filename out of self.spectra and load data using pysynphot.
        Return wave and flux in microns and mJy.
        """
        if self.key not in self.spectra:
            msg = "Provided SED key, %s, not supported." % self.key
            raise EngineInputError(value=msg)

        filepath = os.path.join(default_refdata_directory,
                                "sed",
                                self.__class__.__name__.lower(),
                                self.spectra[self.key]['filename'])
        wave, flux = read_psyn_spectrum(filepath)
        return wave, flux


class HST_CalSpec(Mapped):

    """
    Class implementing catalog of HST calibration stars. Purely mapped with no extra functionality.
    """

    pass


class Brown(Mapped):

    """
    Class implementing catalog of integrated galaxy spectra from Brown et al. (2014).
    If webapp=True, load select subset.
    If webapp=False, load the entire catalog.
    """

    def __init__(self, webapp=False, config={}, **kwargs):
        """
        How an SED is configured depends on the defaults and any input configuration.

        Parameters
        ----------
        webapp: bool
            Determines whether strict API checks should be performed
        config: dict
            Extra configuration information in engine input API dict format
        **kwargs: keyword/value pairs
            Extra configuration information in kwargs format
        """
        SED.__init__(self, webapp=webapp, config=config, **kwargs)
        spectra_file = os.path.join(default_refdata_directory, "sed", self.spectra["config_all"])
        self.spectra = read_json(spectra_file)
        self.wave, self.flux = self.get_spectrum()
        self.wmin = self.wave.min()
        self.wmax = self.wave.max()


def SEDFactory(webapp=False, config={}, **kwargs):
    """
    Function to take configuration data and build/return an appropriately configured instance
    of the desired SED class.

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
    types = recursive_subclasses(SED)
    sed_types = [t.__name__.lower() for t in types]
    type_map = dict(list(zip(sed_types, types)))

    # pop the sed_format out of the input config and error out if it's not there
    try:
        sed_type = all_config["sed_type"]
    except KeyError as e:
        msg = "No SED format specified. (%s)" % e
        raise EngineInputError(value=msg)

    # make sure sed_format is a supported one and build/return appropriate instance if so
    if sed_type not in sed_types:
        msg = "Unsupported or not yet implemented SED: %s" % sed_type
        raise EngineInputError(value=msg)
    else:
        cls = type_map[sed_type](webapp=webapp, config=all_config)
        return cls
