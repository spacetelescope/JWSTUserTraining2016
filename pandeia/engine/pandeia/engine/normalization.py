# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import os
import numpy as np
import astropy.units as u

import pysynphot as psyn

from .config import DefaultConfig
from .custom_exceptions import DataConfigurationError, EngineInputError
from .pandeia_warnings import normalization_warning_messages as warning_messages
from .utils import get_key_list, get_dict_from_keys, recursive_subclasses, merge_data
from .instrument_factory import InstrumentFactory
from .constants import pandeia_waveunits, pandeia_fluxunits

from . import io_utils as io
from . import config as cf

default_refdata_directory = cf.default_refdata_directory


# we use the plural of 'microns' everywhere, but astropy doesn't recognize it. configure it here...
MICRONS = u.def_unit("microns", u.um, format={'generic': 'microns', 'console': 'microns'})
u.add_enabled_units([MICRONS])


class Normalization(DefaultConfig):

    """
    Class for handling configuration and application of spectrum normalization

    Attributes
    ----------
    type: string
        Type of normalization to perform
    norm_fluxunit: string
        Unit for fluxes used in normalization
    norm_waveunit: string
        Unit for wavelength used for reference wavelength, 'norm_wave'
    norm_wave: float
        Reference wavelength in 'norm_waveunit' at which the spectrum will be scaled for type="at_lambda"
    norm_flux: float
        Reference flux at 'norm_wave' to which spectrum is scaled to match for type="at_lambda"
    bandpass: string
        Specifies how bandpass information is looked up or calculated for normalization
        types 'hst', 'jwst', and 'photsys'.
    """

    def __init__(self, webapp=False, config={}, **kwargs):
        """
        How a normalization is configured depends on the defaults, the per-type defaults, and any input parameters.

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

        if not hasattr(self, "type"):
            self.type = self.__class__.__name__.lower().replace('normalize', '')

        if self.type in self.types:
            type_conf = self.types[self.type]
            if "defaults" in type_conf:
                type_defaults = type_conf.pop("defaults")
            else:
                type_defaults = {}
            type_conf = merge_data(type_conf, type_defaults, config, dict(**kwargs))
            self.__dict__ = merge_data(self.__dict__, type_conf)

        # do some sanity checks to make sure inputs make sense
        try:
            self._sanity_checks()
        except AttributeError as e:
            self.warnings['no_sanity_check'] = warning_messages['no_sanity_check'] % (self.__class__.__name__, e)

        # we need to run the API checks here after merging in the normalization type-specific defaults
        if webapp:
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
        ref_dir = os.path.join(default_refdata_directory, "normalization")
        config = io.read_json(os.path.join(ref_dir, "config.json"), raise_except=True)

        # pop the defaults entry out
        if "defaults" in config:
            defaults = config.pop("defaults")
            config.update(defaults)
        else:
            msg = "No normalization defaults defined."
            raise DataConfigurationError(value=msg)

        return config

    def _sanity_checks(self):
        """
        Make sure norm_fluxunit is supported and type is a configured normalization type
        """
        if self.norm_fluxunit not in self.fluxunits:
            msg = "Unsupported flux unit %s for normalization type %s" % (self.norm_fluxunit, self.type)
            raise EngineInputError(value=msg)

        if self.type not in self.types:
            msg = "No configuration data found for normalization type %s" % self.type
            raise DataConfigurationError(value=msg)

    def to_pysynphot(self, wave, flux):
        """
        Take wavelength and flux arrays from pandeia and pack them into a pysynphot Spectrum object.
        Then convert it to self.norm_fluxunit.

        Parameters
        ----------
        wave: 1D np.ndarray
            Vector of wavelengths in pandeia_waveunits
        flux: 1D np.ndarray
            Vector of fluxes in pandeia_fluxunits

        Returns
        -------
        sp: pysynphot Spectrum object
            Spectrum in pysynphot format
        """
        sp = psyn.ArraySpectrum(wave, flux, waveunits=pandeia_waveunits, fluxunits=pandeia_fluxunits)
        sp.convert(self.norm_fluxunit)
        return sp

    def from_pysynphot(self, sp):
        """
        Take pysynphot Spectrum, convert it to pandeia_waveunits and pandeia_fluxunits for pandeia use,
        and return the wavelength and flux vectors.

        Parameters
        ----------
        sp: pysynphot Spectrum object
            Spectrum in pysynphot format

        Returns
        -------
        wave, flux:  1D np.ndarray, 1D np.ndarray
            Wavelength and flux vectors in 'pandeia_waveunits' and 'pandeia_fluxunits', respectively
        """
        sp.convert(pandeia_waveunits)
        sp.convert(pandeia_fluxunits)
        wave, flux = sp.wave, sp.flux
        return wave, flux


class NormalizeAtLambda(Normalization):

    """
    Subclass for handling normalization via a reference flux, self.norm_flux, at a reference wavelength, self.norm_wave.
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
        Normalization.__init__(self, webapp=webapp, config=config, **kwargs)

        # check sanity of input units
        if self.norm_waveunit not in self.waveunits:
            msg = "Unsupported reference wavelength unit, %s.  Supported units are: %s" % (self.norm_waveunit, repr(self.waveunits))
            raise EngineInputError(value=msg)

        if self.norm_fluxunit not in self.fluxunits:
            msg = "Unsupported flux unit, %s. Supported units are: %s" % (self.norm_fluxunit, repr(self.fluxunits))
            raise EngineInputError(value=msg)

        # use astropy.units to convert norm_wave to pandeia_waveunits
        try:
            self.norm_wave = (self.norm_wave * u.Unit(self.norm_waveunit)).to(u.Unit(pandeia_waveunits)).value
        except Exception as e:
            msg = "Failed to convert reference wavelength units, %s, to pandeia wavelength units, %s: %s" % \
                (self.norm_waveunit, pandeia_waveunits, e)
            raise EngineInputError(value=msg)
        self.norm_waveunit = pandeia_waveunits

        # use pysynphot to convert self.norm_flux to pandeia_fluxunits, if necessary.
        if self.norm_fluxunit != pandeia_fluxunits:
            sp = psyn.ArraySpectrum(
                np.array([self.norm_wave]),
                np.array([self.norm_flux]),
                fluxunits=self.norm_fluxunit,
                waveunits=pandeia_waveunits
            )
            sp.convert(pandeia_fluxunits)
            self.norm_flux = sp.flux[0]
            self.norm_fluxunit = pandeia_fluxunits

    def normalize(self, wave, flux):
        """
        Normalize spectrum by scaling it at self.norm_wave.

        Parameters
        ----------
        wave: 1D np.ndarray
            Wavelength vector in 'pandeia_waveunits'
        flux: 1D np.ndarray
            Flux vector in 'pandeia_fluxunits' containing the spectrum

        Returns
        -------
        scaled_wave, scaled_flux: 1D np.ndarray
            Wavelength and scaled flux vectors in pandeia_waveunits and pandeia_fluxunits
        """
        if self.norm_wave > wave.max() or self.norm_wave < wave.min():
            msg = "Specified normalization wavelength, %.2f, not within wavelength bounds of spectrum: (%.2f, %.2f)" % \
                (self.norm_wave, wave.min(), wave.max())
            raise EngineInputError(value=msg)

        spec_flux = np.interp(self.norm_wave, wave, flux)
        if spec_flux == 0.0:
            key = 'normalized_to_zero_flux'
            self.warnings[key] = warning_messages[key]
            scale_factor = 1.0
        else:
            scale_factor = self.norm_flux / spec_flux

        scaled_wave = np.copy(wave)
        scaled_flux = flux * scale_factor
        return scaled_wave, scaled_flux


class NormalizePhotsys(Normalization):

    """
    Subclass for normalizing a spectrum to a reference flux/magnitude in a filter in a supported photometric system
    (e.g. Cousins, Bessel, SDSS)
    """

    def _get_bandpass(self, *args):
        """
        Parse a self.bandpass of the form <photsys>,<filter> and use the keys to look up a bandpass file
        in self.photsystems.  Use pysynphot to load the bandpass and return it.

        Returns
        -------
        bp: pysynphot.SpectralElement
            Bandpass throughput data as returned by pysynphot.FileBandpass() and converted to pandeia
            wavelength units
        """
        keys = get_key_list(self.bandpass, separator=',')  # should be [photsys, filter]
        bp_filename = get_dict_from_keys(self.bandpasses, keys)['filename']
        bp_path = os.path.join(default_refdata_directory, "normalization", "bandpass", bp_filename)
        bp = psyn.FileBandpass(bp_path)
        bp.convert(pandeia_waveunits)
        return bp

    def normalize(self, wave, flux):
        """
        Normalize a spectrum to a bandpass fetched from self._get_bandpass()

        Parameters
        ----------
        wave: 1D np.ndarray
            Wavelength vector in pandeia wavelength units
        flux: 1D np.ndarray
            Flux vector in 'self.norm_fluxunit' containing the spectrum

        Returns
        -------
        scaled_wave, scaled_flux: 1D np.ndarrays
            Wavelength and scaled flux vectors in pandeia units, microns and mJy
        """
        sp = self.to_pysynphot(wave, flux)
        bp = self._get_bandpass(wave)
        # TEMPORARY SHIM #
        # we need to specifically convert the spectrum and bandpass to angstroms
        # to work around a pysynphot unit handling bug
        sp.convert('angstroms')
        bp.convert('angstroms')
        sp_rn = sp.renorm(self.norm_flux, self.norm_fluxunit, bp)
        scaled_wave, scaled_flux = self.from_pysynphot(sp_rn)
        return scaled_wave, scaled_flux


class NormalizeObsmode(NormalizePhotsys):

    """
    Subclass for normalizing a spectrum based on a pysynphot-compatible obsmode string via pysynphot.ObsBandpass()
    """

    def _get_bandpass(self, *args):
        """
        Wrap pysynphot.ObsBandpass() to generate a bandpass based on a valid obsmode specification

        Returns
        -------
        bp: pysynphot.SpectralElement
            pysynphot bandpass converted to pandeia wavelength units.
        """
        if self.webapp and self.bandpass not in self.bandpasses:
            key = "unsupported_normalization_bandpass"
            self.warnings[key] = warning_messages[key] % self.bandpass

        bp = psyn.ObsBandpass(str(self.bandpass))
        bp.convert(pandeia_waveunits)
        return bp


class NormalizeHst(NormalizeObsmode):

    """
    Largely an alias for NormalizeObsmode since HST normalization bands will be specified as obsmode strings.
    """

    pass


class NormalizeJwst(NormalizePhotsys):

    """
    Subclass for normalizing a spectrum based on a JWST configuration
    """

    def _get_bandpass(self, wave):
        """
        Get JWST instrument, mode, filter from self.bandpass, use that to create a JWSTInstrument instance, and then
        query that for throughput vs wavelength.  Use pysynphot.ArrayBandpass() to convert the results to
        pysynphot-compatible form.

        Parameters
        ----------
        wave: 1D np.ndarray
            Wavelength vector onto which JWST bandpass is interpolated

        Returns
        -------
        bp: pysynphot.SpectralElement
            pysynphot bandpass converted to pandeia wavelength units.
        """
        keys = get_key_list(self.bandpass, separator=',')
        if len(keys) != 3:
            msg = "JWST bandpass specification must be of the form <instrument>,<mode>,<filter>"
            raise EngineInputError(value=msg)

        instrument, mode, filt = keys
        config = {}
        config['instrument'] = {}
        config['instrument']['instrument'] = instrument
        config['instrument']['mode'] = mode
        config['instrument']['filter'] = filt
        inst = InstrumentFactory(config=config)
        thruput = inst.get_total_eff(wave)
        bp = psyn.ArrayBandpass(wave, thruput, waveunits=pandeia_waveunits)
        return bp


class NormalizeNone(Normalization):

    """
    Subclass to make no normalization look like the other methods.
    """

    def normalize(self, wave, flux):
        """
        Implement this method to maintain consistent API and simply return what's given.

        Parameters
        ----------
        wave: 1D np.ndarray
            Wavelength vector in pandeia wavelength units
        flux: 1D np.ndarray
            Flux vector in 'self.norm_fluxunit' containing the spectrum

        Returns
        -------
        wave, flux: 1D np.ndarrays
            Return wave and flux without modification
        """
        return wave, flux

    def _sanity_checks(self):
        """
        'type' is the only parameter this normalization type accepts. Should never get here
        if NormalizationFactory is used...
        """
        if self.type not in self.types:
            msg = "No configuration data found for normalization type %s" % self.type
            raise DataConfigurationError(value=msg)


def NormalizationFactory(webapp=False, config={}, **kwargs):

    """
    Function to take normalization configuration data and build/return a configured instance of the appropriate
    subclass of Normalization.

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
    if 'type' in all_config:
        method = all_config['type'].replace('_', '')
    else:
        msg = "Must specify type of normalization to perform."
        raise EngineInputError(value=msg)

    types = recursive_subclasses(Normalization)
    methods = [t.__name__.lower().replace('normalize', '') for t in types]
    type_map = dict(list(zip(methods, types)))

    if method not in methods:
        msg = "Unsupported or not yet implemented normalization method: %s" % method
        raise EngineInputError(value=msg)
    else:
        cls = type_map[method](webapp=webapp, config=config, **kwargs)
        return cls
