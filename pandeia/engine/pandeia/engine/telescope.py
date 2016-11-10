# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import os

from . import io_utils as io
from . import config as cf
from .config import DefaultConfig
from .io_utils import ref_data_interp, ref_data_column
from .pandeia_warnings import telescope_warning_messages as warning_messages

default_refdata_directory = cf.default_refdata_directory


class TelescopeConfig(DefaultConfig):

    """
    This is a trivial top-level class for now. It overrides DefaultConfig._get_config() to get the configuration
    information from a different place, currently under $pandeia_refdata/<tel_name>/telescope.
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
        self.tel_name = self.__class__.__name__.lower()
        self.ref_dir = os.path.join(default_refdata_directory, self.tel_name, "telescope")
        config = io.read_json(os.path.join(self.ref_dir, "config.json"), raise_except=True)
        # add separate CR config, if it's there...
        cr_config = io.read_json(os.path.join(self.ref_dir, "cr_config.json"), raise_except=False)
        config.update(cr_config)
        return config


class Telescope(TelescopeConfig):

    """
    The telescope class contains information about the optical system and other effects
    that are common to all associated instruments.
    """

    def get_ote_eff(self, wave):
        """
        Get efficiency of the telescope optics, the OTE

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate OTE efficiency onto

        Returns
        -------
        ote: numpy.ndarray or float
            If file exists, return ote_efficiency(wave). Else return 1.0.
        """
        ote_file = os.path.join(self.ref_dir, self.paths['otefile'])
        if ote_file is not None:
            ote = ref_data_interp(ote_file, wave, colname='throughput')
        else:
            ote = 1.0
            key = "telescope_ote_efficiency_missing"
            self.warnings[key] = warning_messages[key]
        return ote

    def get_bg_tel(self):
        """
        Get the telescope-specific background

        Parameters
        ----------
        wave: numpy.ndarray
            Wavelength vector to interpolate the telescope background onto

        Returns
        -------
        bg_tel: numpy.ndarray or float
            If file exists, return bg_tel(wave). Else return 0.0
        """
        bg_tel_file = os.path.join(self.ref_dir, self.paths['bg_tel'])
        if bg_tel_file is not None:
            bg_tel_wave = ref_data_column(bg_tel_file, colname='WAVELENGTH')
            bg_tel = ref_data_column(bg_tel_file, colname='SB')
        else:
            bg_tel_wave = [0.1, 30.0]
            bg_tel = [0.0, 0.0]
            key = "telescope_background_missing"
            self.warnings[key] = warning_messages[key]
        return bg_tel_wave, bg_tel
