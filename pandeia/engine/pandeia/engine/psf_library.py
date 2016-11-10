# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import os
import copy

import astropy.io.fits as fits
import numpy as np

from .custom_exceptions import DataError


class PSFLibrary:

    """
    Class for encapsulating a hard PSF library, typically one generated
    by WebbPSF or its derivatives.

    Parameters
    ----------
    path : string, optional
        the path to the PSF library root (any .fits files in the
        path root will be assumed to be a valid PSF file in the WebbPSF format.
        If it is none, the PSF library from the package distribution will be used.
    aperture : string, optional
        restrict psfs to a specific mode with the name matching the keyword string.
        Strictly, the reader will look for psf file names with a substring matching the aperture
        string. 'all' reads all psfs in the path. This is the default. The aperture string in the filename
        should match the .fits keyword APERTURE.

    Methods
    -------
    read_library
    get_values
    get_good_psfs
    get_psf

    """

    def __init__(self, path=None, aperture='all'):
        if path is None:
            path = os.path.join(os.path.dirname(__file__), 'refdata', 'psfs')
        self.read_library(path, aperture)

    def read_library(self, path, aperture, wave_unit='m'):
        """
        Read in a library of PSFs stored in FITS files. The filenames contain information about each PSF. The fields are
        separated by _'s with the resulting list entries corresponding to:
            0 - instrument name
            1 - apertures for which PSF is used. this is a - separated string that is split() into a list of apertures
            2 - wavelength at which PSF is calculated
            3, 4 (optional) - r,theta polar coordinates in FOV at which at which PSF is calculated

        Parameters
        ----------
        path: string
            Path containing the PSF FITS images

        wave_unit: string
            Wavelength unit used by PSF library (default: m)
        """
        files = os.listdir(path)

        psf_files = []
        for file in files:
            filename, ext = os.path.splitext(file)
            filename_split = filename.split('_')
            filename_instrument = filename_split[0]
            filename_apertures = filename_split[1].split('-')
            filename_wave = filename_split[2]
            # If the psf is on-axis, the offset polar coordinates may not exist
            try:
                filename_radius = filename_split[3]
                filename_theta = filename_split[4]
            except:
                filename_radius = 0.
                filename_theta = 0.

            if aperture is 'all':
                if ext == '.fits':
                    psf_files.append(file)
            else:
                if ext == '.fits' and aperture in filename_apertures:
                    psf_files.append(file)

        self._psfs = []

        for psf_file in psf_files:
            hdulist = fits.open(path + '/' + psf_file, memmap=False)

            psf_int = copy.deepcopy(hdulist[0].data)

            ins = hdulist[0].header['INSTRUME'].lower()
            nwaves = hdulist[0].header['NWAVES']
            if nwaves != 1:
                raise ValueError("Input PSF must be monochromatic (nwaves = #r)" % nwaves)

            wave_scl = self._get_unit_scl(wave_unit)
            wave = hdulist[0].header['WAVE0'] * wave_scl

            pix_scl = hdulist[0].header['PIXELSCL']
            diff_limit = hdulist[0].header['DIFFLMT']
            aperture_name = hdulist[0].header['APERTURE'].lower()
            # offset_r/offset_t is the offset from the optical center in webbpsf.  however, for NIRCam masklwb
            # and maskswb, this includes the offset from the optical center to the optimal position along the bar.
            # what we need for source association is instead the offset from that optimal position. this is stored in the
            # optoff_r/optoff_t keywords which need to be provided in these cases. aperture_name here will have these
            # aperture and filter concatenated so do a string compare to look for aperture.
            if 'masklwb' in aperture_name or 'maskswb' in aperture_name:
                try:
                    radius = hdulist[0].header['OPTOFF_R']
                    theta = hdulist[0].header['OPTOFF_T']
                except KeyError as e:
                    msg = "PSFs for bar-shaped masks require offsets from optimal bar position. (%s)" % e
                    raise DataError(value=msg)
            else:
                radius = hdulist[0].header['OFFSET_R']
                theta = hdulist[0].header['OFFSET_T']

            upsamp = hdulist[0].header['DET_SAMP']

            psf = {
                'int': psf_int,
                'wave': wave,
                'pix_scl': pix_scl,
                'diff_limit': diff_limit,
                'upsamp': upsamp,
                'instrument': ins,
                'aperture_name': aperture_name,
                'source_offset': (radius, theta)
            }

            self._psfs.append(psf)
            hdulist.close()

    def get_values(self, key, instrument, aperture_name, source_offset=(0, 0)):
        """
        Returns the available values of a given key for a given instrument and aperture name.

        Parameters
        ----------
        key: string
            Desired key. One of 'int', 'wave', 'pix_scl', 'diff_limit', 'upsamp', 'instrument',
            'aperture_name', 'source_offset_r', or 'source_offset_theta'

        Returns
        -------
        result: tuple of format (list, list)
            Lists of index id's and values
        """
        values = [psf[key] for psf in self._psfs
                  if (instrument == psf['instrument'] and
                      aperture_name in psf['aperture_name']) and
                      source_offset == psf['source_offset']]
        ids = [i for i, psf in enumerate(self._psfs)
               if (instrument == psf['instrument'] and
                   aperture_name in psf['aperture_name']) and
                   source_offset == psf['source_offset']]
        result = ids, values
        return result

    def get_psf(self, wave, instrument, aperture_name, source_offset=(0, 0)):
        """
        Get PSF given wavelength, instrument, instrument aperture, and optionally source offset polar coordinates.
        Values are interpolated linearly between the two PSF library entries that bracket the given wavelength.

        Parameters
        ----------
        wave: float
            Desired wavelength for the PSF
        instrument: str
            Instrument name
        aperture_name: str
            Name of the instrument aperture

        Returns
        -------
        psf: dict
            Dict containing the PSF and associated information
        """
        wids, psf_waves = self.get_values('wave', instrument, aperture_name, source_offset=source_offset)
        nids, nearest_waves = self._find_two_nearest(psf_waves, wave)
        ids = [wids[nids[0]], wids[nids[1]]]

        pix_scl0 = self._psfs[ids[0]]['pix_scl']
        pix_scl1 = self._psfs[ids[1]]['pix_scl']
        if pix_scl0 != pix_scl1:
            raise ValueError("Pixel scales in the library must be the same for a single instrument aperture.")

        psf_int0 = self._psfs[ids[0]]['int']
        psf_int1 = self._psfs[ids[1]]['int']
        wave0 = self._psfs[ids[0]]['wave']
        wave1 = self._psfs[ids[1]]['wave']
        diff_limit0 = self._psfs[ids[0]]['diff_limit']
        diff_limit1 = self._psfs[ids[1]]['diff_limit']

        upsamp0 = self._psfs[ids[0]]['upsamp']
        upsamp1 = self._psfs[ids[1]]['upsamp']
        if upsamp0 != upsamp1:
            raise ValueError("Upsampling factors in the library must be the same for a single instrument aperture.")

        psf_int = psf_int0 + (psf_int1 - psf_int0) * (wave - wave0) / (wave1 - wave0)
        diff_limit = diff_limit0 + (diff_limit1 - diff_limit0) * (wave - wave0) / (wave1 - wave0)

        psf = {
            'int': psf_int,
            'wave': wave,
            'pix_scl': pix_scl0,
            'diff_limit': diff_limit,
            'upsamp': upsamp0,
            'instrument': self._psfs[ids[0]]['instrument'],
            'aperture_name': self._psfs[ids[0]]['aperture_name'],
            'source_offset': self._psfs[ids[0]]['source_offset']
        }
        return psf

    def get_pix_scale(self, instrument, aperture_name):
        """
        Get PSF pixel scale for given instrument/aperture

        Parameters
        ----------
        instrument: str
            Instrument name
        aperture_name: str
            Name of instrument aperture

        Returns
        -------
        pix_scl: float
            Pixel scale of the PSF in arcsec/pixel
        """
        wids, psf_waves = self.get_values('wave', instrument, aperture_name)
        pix_scl = self._psfs[wids[0]]['pix_scl']
        return pix_scl

    def get_shape(self, instrument, aperture_name):
        """
        Get PSF shape for given instrument/aperture

        Parameters
        ----------
        instrument: str
            Instrument name
        aperture_name: str
            Name of instrument aperture

        Returns
        -------
        sh: tuple (int, int)
            Shape of the PSF kernel image
        """
        wids, psf_waves = self.get_values('wave', instrument, aperture_name)
        sh = self._psfs[wids[0]]['int'].shape
        return sh

    def get_upsamp(self, instrument, aperture_name):
        """
        Get PSF upsampling for given instrument/aperture

        Parameters
        ----------
        instrument: str
            Instrument name
        aperture_name: str
            Name of instrument aperture

        Returns
        -------
        upsamp: int
            PSF upsampling factor
        """
        wids, psf_waves = self.get_values('wave', instrument, aperture_name)
        upsamp = self._psfs[wids[0]]['upsamp']
        return upsamp

    def get_offsets(self, instrument, aperture_name):
        """
        Get available PSF offset positions. PSF libraries supporting position-dependent
        PSFs will have multiple. Libraries that do not support position-dependent PSFs will
        only have one position.

        Parameters
        ----------
        list: list-like
            List of values to search
        value: float
            Input value to be bracketed

        Returns
        -------
        psf_positions: list of tuples (source_offset_r, source_offset_theta)
            Polar coordinates of PSF offsets.

        """
        offsets = [psf['source_offset'] for psf in self._psfs if instrument ==
                   psf['instrument'] and aperture_name in psf['aperture_name']]
        unique_offsets = list(set(offsets))
        return unique_offsets

    def associate_offset_to_source(self, sources, instrument, aperture_name):

        psf_offsets = self.get_offsets(instrument, aperture_name)
        psf_associations = []
        for source in sources:
            source_offset_radius = np.sqrt(source.position['x_offset']**2. + source.position['y_offset']**2.)
            # Currently, we only associate radius, not angle.
            distances = [np.abs(source_offset_radius-psf_offset[0]) for psf_offset in psf_offsets]
            closest_index = np.argmin(distances)
            psf_associations.append(psf_offsets[closest_index])

        return psf_associations

    def _find_two_nearest(self, list, value):
        """
        Find the subscripts of the two neighboring values in list relative to
        an input value.

        Parameters
        ----------
        list: list-like
            List of values to search
        value: float
            Input value to be bracketed

        Returns
        -------
        vals: tuple (list, np.ndarray)
            Bracketing indices in both list and numpy array format
        """
        nplist = np.array(list)

        # ensure that value is not outside the range of list
        if value < nplist.min() or value > nplist.max():
            raise ValueError("Value must be within range of list.")

        diff = np.abs(nplist - value)
        diff_sorted = np.sort(diff)

        close = diff_sorted[1]
        closer = diff_sorted[0]

        ids = [i for i, v in enumerate(diff) if diff[i] == close or diff[i] == closer]

        vals = ids, nplist[ids]
        return vals

    def _get_unit_scl(self, unit):
        """
        Get unit scale factor for given unit string: m, micron, nm, or Angstrom.
        This should get replaced with real astropy.units support eventually.

        Parameters
        ----------
        unit: str
            Unit string

        Returns
        -------
        scale: float
            Scale factor corresponding to unit string
        """
        scales = {'m': 1.0e6, 'micron': 1.0, 'nm': 1.0e-3, 'Angstrom': 1.0e-4}
        try:
            scale = scales[unit]
            return scale
        except KeyError:
            raise KeyError('Unknown wavelength unit')
