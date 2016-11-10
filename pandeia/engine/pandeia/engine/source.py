# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import os

from .config import DefaultConfig
from .custom_exceptions import DataConfigurationError, EngineInputError
from .utils import merge_data

from . import config as cf
from . import io_utils as io

default_refdata_directory = cf.default_refdata_directory


class Line(DefaultConfig):

    """
    Container for spectral line parameters

    Parameters
    ----------
    config: dict (optional)
        dictionary containing necessary configuration information
    **kwargs: list of keyword/value pairs
        parameter keyword and value pairs to augment defaults and config

    Attributes
    ----------
    id: str or int
        Line ID tag
    name: str
        Name given to the line
    center: float
        Central wavelength of the line in microns
    width: float
        Line FWHM in km/s
    strength: float
        Line strength in erg/cm^2/s for emission or central optical depth for absorption
    profile: str
        Line profile change, currently only "gaussian" is supported
    emission_or_absorption: str
        Line profile type, either 'emission' or 'absorption'
    """
    def _get_config(self):
        """
        Read default configuration from JSON

        Returns
        -------
        config: dict
            All desired class attributes should have defaults defined in the config file
        """
        # use this trick to key the configuration file name off the name of the instantiated subclass
        ref_dir = os.path.join(default_refdata_directory, "source", "line")
        config = io.read_json(os.path.join(ref_dir, "defaults.json"), raise_except=True)
        return config


class Source(DefaultConfig):

    """
    Container for source parameters

    Parameters
    ----------
    config: dict (optional)
        dictionary containing necessary configuration information
    **kwargs: list of keyword/value pairs
        parameter keyword and value pairs to augment defaults and config

    Attributes
    ----------
    id: str or int
        source id
    position: dict
        x_offset: float
            x offset (arcsec)
        y_offset: float
            y offset (arcsec)
        orientation: float
            orientation angle (degrees)
    shape: dict
        geometry: string
            supported source geometries are:
                "point" - requires no further parameters

                "flat" - requires the following parameters:
                    major: float
                        major axis length (arcsec).
                    minor: float
                        minor axis length (arcsec).

                "gaussian2d" - requires "major" and "minor"

                "sersic" - requires "major", "minor", and an additional parameter:
                    sersic_index: float
                        power law index that sets the shape of a sersic profile. used only for "sersic" geometry.
                            sersic_index = 1.0 --> exponential
                            sersic_index = 0.5 --> gaussian
                            sersic_index = 4.0 --> de Vaucouleurs
    spectrum: dict
        lines: list
            list of line definitions
        name: string
            Name of the source
        redshift: float
            Redshift of the source
        sed: dict
          Defines the spectral energy distribution of the spectrum.

            sed_type: string
                Type of the spectral energy distribution. Each type requires its own set
                of parameters. The analytic sed_type's (none, flat, powerlaw, flat) all
                require 'wmin', 'wmax', and 'sampling' to define the range and wavelength
                sampling over which the model spectrum is calculated. However, they are only
                available in the API for testing purposes and should not be configured via
                the UI.

                    no_continuum - No continuum, specifically F(lam) = 0.0 over specified range [wmin, wmax]
                        wmin: float (default 0.5)
                            Minimum wavelength in microns
                        wmax: float (default 30.0)
                            Maximum wavelength in microns
                        sampling: int (default 200)
                            Sets the logarithmic wavelength sampling of the model spectrum

                    flat - Flat spectrum in specified units calculated over specified range [wmin, wmax]
                        wmin: float (default 0.5)
                            Minimum wavelength in microns
                        wmax: float (default 30.0)
                            Maximum wavelength in microns
                        sampling: int (default 200)
                            Sets the logarithmic wavelength sampling of the model spectrum
                        unit: string
                            Units of spectrum, either 'fnu' or 'flam'

                    powerlaw - Powerlaw spectrum where F ~ lambda ^ index calculated over range [wmin, wmax]
                        wmin: float (default 0.5)
                            Minimum wavelength in microns
                        wmax: float (default 30.0)
                            Maximum wavelength in microns
                        sampling: int (default 200)
                            Sets the logarithmic wavelength sampling of the model spectrum
                        unit: string
                            Units of spectrum, either 'fnu' or 'flam'
                        index: float
                            Exponent of the power law

                    blackbody - Blackbody spectrym calculated over range [wmin, wmax]
                        wmin: float (default 0.5)
                            Minimum wavelength in microns
                        wmax: float (default 30.0)
                            Maximum wavelength in microns
                        sampling: int (default 200)
                            Sets the logarithmic wavelength sampling of the model spectrum
                        temp: float
                            Temperature of the blackbody

                    phoenix - Parameterized stellar atmosphere models calculated by the Phoenix group
                        key: string
                            In webapp mode, a key is used to look up a predefined set of parameters. If not
                            in webapp mode and if key is not provided, model parameters can be passed directly:
                        teff: float
                            Effective temperature. Allowed range is 2000 K to 70000 K
                        log_g: float
                            Logarithm of the surface gravity in cgs units. Allowed range is 0.0 to 5.5.
                        metallicity: float
                            Logarithm of the metallicity in units of solar metallicity. Allowed range is -4.0 to 0.5.

                    hst_calspec - HST standard star spectra
                        key: string
                            Key used to look up which spectrum load.

                    galaxies - Integrated spectra of galaxies from Brown et al. (2014)
                        key: string
                            Key used to look up which spectrum load.

                    input - spectrum provided via input arrays
                        spectrum: list-like or numpy.ndarray
                            The 0th index is taken to be wavelength in units of 'mJy'.
                            The 1st index is taken to be the flux in units of 'microns'.

        normalization: dict
            Define the brightness of the source
                type: string
                    Method of normalization to perform.
                    Supported methods are: "at_lambda", "hst", "jwst", "photsys", and "none"
                norm_wave: float
                    Reference wavelength in 'norm_waveunit' at which spectrum will be scaled for type 'at_lambda'.
                norm_flux: float
                    Reference flux in 'norm_fluxunit' to which spectrum will be scaled.
                norm_fluxunit: string
                    Specify the flux units in which the normalization should occur.
                    Supports flam, fnu, vegamag, abmag, mjy, ujy, njy, jy
                norm_waveunit: string
                    Specify the wavelength units used in normalization
                bandpass: string
                    Specifies how bandpass information is looked up or calculated for normalization
                    types 'hst', 'jwst', and 'photsys'.

    Methods
    -------
    add_line, remove_line
    """
    def __init__(self, config={}, webapp=False, **kwargs):
        """
        How a source is configured depends on the defaults, the shape-specific defaults, and any input parameters.

        Parameters
        ----------
        webapp: bool
            Determines whether strict API checks should be performed
        config: dict
            Extra configuration information in engine input API dict format
        **kwargs: keyword/value pairs
            Extra configuration information in kwargs format
        """
        DefaultConfig.__init__(self, webapp=webapp, config=config, **kwargs)

        if self.shape['geometry'] in self.geometries:
            shape_config = self.geometries[self.shape['geometry']]
            # clean out stuff not needed anymore by engine
            shape_config.pop("display_string")
            self.shape = merge_data(shape_config, self.shape)
        else:
            msg = "Unsupported source geometry: %s" % self.shape['geometry']
            raise EngineInputError(value=msg)

        all_config = merge_data(config, dict(**kwargs))

        if webapp:
            self._nested_api_checks(all_config)

    def _get_config(self):
        """
        Read default configuration from JSON

        Returns
        -------
        config: dict
            All desired class attributes should have defaults defined in the config file
        """
        # use this trick to key the configuration file name off the name of the instantiated subclass
        ref_dir = os.path.join(default_refdata_directory, "source")
        config = io.read_json(os.path.join(ref_dir, "config.json"), raise_except=True)

        return config

    def _nested_api_checks(self, conf):
        """
        If webapp=True, check to make sure all of the configured shape_parameters and spectrum_parameters exist in the input
        configuration.  api_parameters only works at top level, but we need to check one level down
        in self.shape and self.spectrum which is why this is a separate check.  Proper fix would be to refactor Source
        into sub-classes for spectrum and each geometry type.
        """
        # 'shape' and 'spectrum' are already checked for via api_parameters
        for section in ['shape', 'spectrum', 'position']:
            if section in conf:
                parskey = '%s_parameters' % section
                params = getattr(self, section)[parskey]
                if parskey in getattr(self, section):
                    for p in params:
                        if p not in conf[section]:
                            key = "source_%s_parameter_%s_missing" % (section, p)
                            msg = "Source %s configuration missing API parameter, %s. Using the default value of %s." % (
                                section,
                                p,
                                getattr(self, section)[p]
                            )
                            self.warnings[key] = msg
                else:
                    key = "source_%s_parameters_missing" % section
                    msg = "Missing list of required API parameters for source %s." % section
                    self.warnings[key] = msg

    def _sanity_checks(self):
        """
        Make sure normalization, SED, and shape information are defined.
        """
        # need this information to make a source...
        if 'normalization' not in self.spectrum:
            msg = "Normalization information missing for default source."
            raise DataConfigurationError(value=msg)
        if 'sed' not in self.spectrum:
            msg = "SED information missing for default source."
            raise DataConfigurationError(value=msg)
        if not hasattr(self, "shape"):
            msg = "Shape information missing for default source."
            raise DataConfigurationError(value=msg)
        # make sure configured shape_parameters exist in the data files
        if "shape_parameters" in self.shape:
            for p in self.shape['shape_parameters']:
                if p not in self.shape:
                    msg = "Shape information is missing %s for default source." % p
                    raise DataConfigurationError(value=msg)

    def add_line(self, config={}, **kwargs):
        """
        Add a spectral Line to the Source's spectrum.
        """
        # use a DefaultConfig subclass here for convenience to populate
        # defaults, but store the dict format in spectrum['lines'] to simplfify
        # instantiation from JSON
        line_definition = Line(config=config, **kwargs)
        self.warnings.update(line_definition.warnings)

        self.spectrum['lines'].append(line_definition.__dict__)

    def remove_line(self, line_definition):
        """
        Remove a given line_definition from Source's spectrum
        """
        self.spectrum['lines'].remove(line_definition)
