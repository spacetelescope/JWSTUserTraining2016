# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import numpy as np
import scipy.constants as cs
import scipy.integrate as ig
import pyfftw
from astropy.io import fits
from astropy.convolution import convolve_fft

from .normalization import NormalizationFactory
from .extinction import ExtinctionFactory
from .sed import SEDFactory
from .coords import Grid
from .utils import merge_wavelengths, spectrum_resample
from .custom_exceptions import EngineInputError, WavesetMismatch, DataError, RangeError
from .pandeia_warnings import astrospectrum_warning_messages as warning_messages

from six.moves import range

# pyfftw.interfaces.cache.enable()


class ModelSceneCube(object):

    """
    This generates a model intensity cube of a scene. It takes the analytic descriptions
    of the components of the scene and samples them into a Grid.

    Here we put wavelength as the 3rd index to enable broadcasting. The cube is a stack
    of images, each at a different wavelength. If we wish to add or multiply a wavelength
    vector (e.g. throughput) with the cube, the cube needs to be in this index order for
    it to work efficiently. FITS, however, expects an order of (wave, y, x) so a transpose
    will be required to convert to that format. Other formats may have similar requirements.

    Parameters
    ----------
    source_spectra: list of pandeia.engine.astro_spectrum.AstroSpectrum instances
        Source spectra to add to the model cube
    grid: pandeia.engine.coord.Grid instance
        Grid describing the 2D spatial coordinates for each plane of the cube

    Methods
    -------
    add_source:  add an AstroSpectrum to the cube
    export_to_fits: export model cube to a FITS file
    _sersic_profile:  generate a source intensity distribution based on a Sersic intensity profile
    _point_source: generate a point source with sub-pixel positioning
    """

    def __init__(self, source_spectra, grid):
        self.grid = grid
        self.x = grid.x
        self.y = grid.y
        self.xsamp = grid.xsamp
        self.ysamp = grid.ysamp
        self.nx = grid.nx
        self.ny = grid.ny
        self.wave = source_spectra[0].wave
        # now go through the source spectra and make sure they've been merged
        # onto the same wavelength set
        for s in source_spectra:
            if not np.array_equal(self.wave, s.wave):
                message = "Model cube input spectra must be sampled at the same wavelengths."
                raise EngineInputError(value=message)
        # stick wavelength as the 3rd index of the cube to enable broadcasting
        self.int = np.zeros(self.x.shape + (self.wave.size,))

        for source_spectrum in source_spectra:
            self.add_source(source_spectrum)

    def add_source(self, spectrum):
        """
        Add a source to the model cube.

        Parameters
        ----------
        spectrum: pandeia.engine.astro_spectrum.AstroSpectrum instance
            Spectrum object containing spatial and spectral information for the source to be added
        """
        src = spectrum.src
        if src.shape['geometry'] == "point":
            plane = self._point_source(
                xoff=src.position['x_offset'],
                yoff=src.position['y_offset']
            )
        elif src.shape['geometry'] == "gaussian2d":
            plane = self._sersic_profile(
                major=src.shape['major'] * np.sqrt(2.0),  # multiply by sqrt(2) to get gaussian sigma
                minor=src.shape['minor'] * np.sqrt(2.0),
                pa=src.position['orientation'],
                xoff=src.position['x_offset'],
                yoff=src.position['y_offset'],
                sersic_index=0.5
            )
        elif src.shape['geometry'] == "sersic":
            plane = self._sersic_profile(
                major=src.shape['major'],
                minor=src.shape['minor'],
                pa=src.position['orientation'],
                xoff=src.position['x_offset'],
                yoff=src.position['y_offset'],
                sersic_index=src.shape['sersic_index']
            )
        elif src.shape['geometry'] == "flat":
            plane = self._flat_source(
                major=src.shape['major'],
                minor=src.shape['minor'],
                pa=src.position['orientation'],
                xoff=src.position['x_offset'],
                yoff=src.position['y_offset']
            )
        else:
            msg = "Unsupported source geometry: %s" % src.shape['geometry']
            raise EngineInputError(value=msg)

        self.int += plane.reshape(plane.shape + (1,)) * spectrum.flux

    def export_to_fits(self, fitsfile='ModelSceneCube.fits'):
        """
        Write model cube to a FITS file
        """
        header = self.grid.wcs_info()
        fits.writeto(fitsfile, np.rollaxis(self.int, 2), header, clobber=True)
        fits.append(fitsfile, self.wave)

    def _sersic_profile(self, major, minor, pa=0.0, xoff=0.0, yoff=0.0, sersic_index=1.0):
        """
        Create a 2-dimensional elliptical source on the current grid. The intensity profile
        is described by a Sersic profile, I(r) = I(0) * exp(-(r/r_scale)**(1/n)), where r_scale
        is the scale length where I(r) = I(0)/e and n is the Sersic index. The ellipticity
        is governed by specifying major and minor axis scale lengths separately.

        Normalization is performed using scipy.integrate.dblquad to integrate the profile in 2D
        from -Inf to +Inf in both axes.  This will account for flux that falls outside of the FOV.
        This works well for sersic_index <= 2.0 or so, but slow convergence for larger values can impact
        performance noticeably.  For example, it takes 10 times longer to normalize a de Vaucouleurs profile with
        sersic_index=4.0 than an exponential disc with sersic_index=1.  Truncating the integration range can help,
        but also significantly impacts accuracy as sersic_index increases.  Ultimate solution may be to implement
        the integration in C somehow.

        Parameters
        ----------
        minor: float
            Minor axis scale length
        major: float
            Major axis scale length
        pa: float
            Position angle in degrees of major axis measured positive in +X direction
        xoff: float
            Offset in X direction
        yoff: float
            Offset in Y direction
        sersic_index: float
            Sersic profile shape parameter. 0.5 => gaussian, 1.0 => exponential, 4.0 => de Vaucouleurs

        Returns
        -------
        g: 2D numpy.ndarray
            2D image containing normalized source intensity
        """
        yrot, xrot = self.grid.shift_rotate(yoff, xoff, pa)
        xrot = np.abs(xrot)
        yrot = np.abs(yrot)
        g = self._sersic_func(yrot, xrot, major, minor, sersic_index)
        # integrate the Sersic profile to get the total flux for normalization, including flux outside the FOV
        integral = ig.dblquad(
            self._sersic_func,
            -np.inf,
            np.inf,
            lambda y: -np.inf,
            lambda y: np.inf,
            args=(major, minor, sersic_index)
        )
        norm = integral[0] / (self.grid.xsamp * self.grid.ysamp)  # convert area in arcsec to area in pixels
        g = g / norm
        return g

    def _sersic_func(self, y, x, major, minor, index):
        """
        Implement Sersic intensity profile in a way that scipy.integrate.dblquad can use since it is passed there
        as a callable. We need to integrate this function to normalize it properly for flux outside the calculation FOV.

        Parameters
        ----------
        y: float or numpy.ndarray
            Y values for evaluating function
        x: float or numpy.ndarray
            X values for evaluating function
        major: float
            Major axis scale length
        minor: float
            Minor axis scale length
        index: float
            Sersic index

        Returns
        -------
        profile: float or numpy.ndarray
            Float or array containing evaluated Sersic profile
        """
        dist = np.sqrt((x / major)**2.0 + (y / minor)**2.0)
        profile = np.exp(-dist**(1.0 / index))
        return profile

    def _point_source(self, xoff=0.0, yoff=0.0):
        """
        Use pandeia.engine.coords.Grid.point_source() to generate a point source with subpixel
        positioning.

        Parameters
        ----------
        xoff: float
            Offset in X direction
        yoff: float
            Offset in Y direction

        Returns
        -------
        pt_src: 2D numpy.ndarray
            2D image containing normalized point source intensity
        """
        pt_src = self.grid.point_source(xoff=xoff, yoff=yoff)
        return pt_src

    def _flat_source(self, major, minor, pa=0.0, xoff=0.0, yoff=0.0):
        """
        Implement a source with a constant surface brightness as an ellipse, i.e. a tilted disc.

        Parameters
        ----------
        minor: float
            Minor axis scale length
        major: float
            Major axis scale length
        pa: float
            Position angle in degrees of major axis measured positive in +X direction
        xoff: float
            Offset in X direction
        yoff: float
            Offset in Y direction

        Returns
        -------
        src: 2D numpy.ndarray
            2D image containing normalized source intensity
        """
        src = self.grid.elliptical_mask(major, minor, pa=pa, xoff=xoff, yoff=yoff)
        total_area = np.pi * major * minor / (self.xsamp * self.ysamp)  # total area of ellipse in pixels
        # mask values are pixel fraction enclosed by ellipse. divide by area to normalize to fraction of total flux.
        src = src / total_area
        return src


class AstroSpectrum(object):

    """
    Class that contains a 1D wavelength grid and methods to calculate a single point spectrum,
    for instance from a point source.

    Parameters
    ----------
    src: pandeia.engine.source.Source
        Source instance containing parameters relevant for the definition of the astro spectrum.

    Attributes
    ----------
    src: pandeia.engine.source.Source instance
        The input source containing parameters relevant for the definition of the astro spectrum.
    line_definition: list of pandeia.engine.source.Line instances
        The lines contained in the spectrum of self.src
    normalization: instance of pandeia.engine.normalization.Normalization subclass
        Normalization information defining the brightness of self.src
    sed: instance of pandeia.engine.sed.SED subclass
        Spectral energy distribution information for self.src (excluding lines)
    wave: 1D np.ndarray
        The wavelength vector of the spectrum
    flux: 1D np.ndarray
        The flux vector of the spectrum
    nw: int
           Number of wavelength points

    Methods
    -------
    export_to_fits_table()
        Write self.wave and self.flux to binary FITS table
    add_spectral_line()
        Add a spectral line to the spectrum
    resample()
        Resample spectrum to a new set of wavelengths via interpolation
    trim()
        Trim spectrum to a specified wavelength range

    Notes
    -----

    """

    def __init__(self, src, webapp=False, **kwargs):

        self.src = src
        self.warnings = {}
        # get the spectrum information from the source
        self.line_definitions = self.src.spectrum['lines']
        self.normalization = NormalizationFactory(config=self.src.spectrum['normalization'], webapp=webapp)
        self.extinction = ExtinctionFactory(config=self.src.spectrum['extinction'], webapp=webapp)

        # if redshift is provided (it should be if config file is up-to-date) use it, else assume it's 0.0
        if 'redshift' in self.src.spectrum:
            z = self.src.spectrum['redshift']
        else:
            z = 0.0

        if z <= -1.0:
            msg = "Specified source redshift must be > -1.0."
            raise EngineInputError(value=msg)

        self.sed = SEDFactory(config=self.src.spectrum['sed'], webapp=webapp, z=z)

        self.wave, self.flux = (1.0 + z) * self.sed.wave, self.sed.flux
        self.nw = self.wave.size

        # apply extinction after redshift
        # the normalization is passed in so the bandpass can be obtained through its helper functions
        self.wave, self.flux = self.extinction.extinction(self.wave, self.flux)

        # apply normalization after extinction
        self.wave, self.flux = self.normalization.normalize(self.wave, self.flux)

        # update warnings...
        self.warnings.update(self.normalization.warnings)
        self.warnings.update(self.sed.warnings)

        # add lines after normalizing continuum since line strengths are currently specified
        # using physical units. future support for specifying line strengths as equivalent widths
        # would allow lines to be added before normalization.
        if len(self.line_definitions) > 0:
            self.add_spectral_lines()

    def export_to_fits_table(self, fitsfile='source_spectrum.fits'):
        c1 = fits.Column(array=self.wave, format='D', name='wavelength')
        c2 = fits.Column(array=self.flux, format='D', name='flux')
        tbhdu = fits.BinTableHDU.from_columns([c1, c2])
        tbhdu.writeto(fitsfile, clobber=True)

    def add_spectral_lines(self):
        """
        This method adds lines to the spectrum by either adding the line flux (for emission lines) or calculating
        an absorption line using the optical depth, F = Fcont * exp(-tau). The sampling of the spectrum is checked for each line
        and new samples added as necessary to assure adequate sampling for each line. The sampling of the underlying spectrum is
        doubled to avoid any Nyquist smoothing effects.
        """
        # make a nyquist sampled version of the underlying wavelength set to mitigate the effects of resampling.
        # pysynphot is used for the resampling so flux is properly preserved in the process.
        halves = self.wave[0:-1] + 0.5 * np.diff(self.wave)
        nywave = np.append(self.wave, halves)
        nywave.sort()
        lines = []

        # loop through line definitions to build waveset properly sampled for all of them
        for line_definition in self.line_definitions:
            spectral_line = SpectralLine(line_definition)
            lines.append(spectral_line)
            # update wavelength set so that it optimally samples the line's wavelength region.
            nywave = spectral_line.line_waveset(nywave)

        # resample underlying spectrum to new wavelength set.
        # this method will update self.wave, self.flux, and self.nw accordingly.
        self.resample(nywave)

        # now apply the lines to the spectrum
        for spectral_line in lines:
            if spectral_line.emission_or_absorption == 'emission':
                self.flux += spectral_line.flux(self.wave)
            elif spectral_line.emission_or_absorption == 'absorption':
                self.flux *= np.exp(-spectral_line.optical_depth(self.wave))
            else:
                raise EngineInputError(value="Invalid spectral line type: %s" % spectral_line.emission_or_absorption)

    def resample(self, wavelengths):
        """
        Re-sample spectrum onto a new set of wavelengths.

        This method modifies self.nw, self.wave, and self.flux.

        Arguments
        =========
        wavelengths: 1D numpy array
            Array of new wavelengths to sample spectrum onto
        """
        # we use np.nan to fill in fluxes for wavelengths out of range of the original spectra.
        # convert them back to 0.0 since the 2D convolution routine doesn't like NaN's.
        self.flux = np.nan_to_num(spectrum_resample(self.flux, self.wave, wavelengths))
        self.wave = wavelengths
        self.nw = self.wave.size

    def trim(self, wmin, wmax):
        """
        Trim spectrum to wavelength range specified by wmin and wmax

        Parameters
        ----------
        wmin: float
            Minimum wavelength (microns)
        wmax: float
            Maximum wavelength (microns)
        """
        valid_subs = np.where((self.wave >= wmin) & (self.wave <= wmax))
        if valid_subs is None:
            msg = "Spectrum for source does not have any overlap with instrument wavelength range."
            raise WavesetMismatch(value=msg)

        trim_wave = []
        # make sure trimmed wavelengths include wmin and wmax
        if self.wave[valid_subs][0] != wmin:
            trim_wave.append(wmin)

        trim_wave.extend(self.wave[valid_subs])

        if self.wave[valid_subs][-1] != wmax:
            trim_wave.append(wmax)

        trim_wave = np.array(trim_wave)

        # this will use pysynphot to properly sample the endpoints
        trim_flux = spectrum_resample(self.flux, self.wave, trim_wave)

        self.wave = trim_wave
        self.flux = trim_flux
        self.nw = self.wave.size


class SpectralLine(object):

    """
    A spectral line - that is the information and methods needed to add a single line to a spectrum.
    The line could also form a component of a composite line (for instance a P Cygni line consisting
    of an emission line AND an absorption line). In this case, one would create two different instances
    of SpectralLine, with relevant profile parameters, and add them both to the spectrum. The class
    includes methods to create a wavelength window with finer sampling than the rest of the spectrum.
    The single line is modeled as a gaussian profile (only option right now, could be more later)
    with a FWHM, a central wavelength and a strength. The line can be either in emission or
    absorption. In the case of an emission line, the strength is an integrated line intensity in
    cgs units, and in the case of an absorption line the strength is a peak optical depth.

    Arguments
    definition - A dictionary containing the line parameters (center, width, strength, profile and emission_or_absorption).
    """

    def __init__(self, definition):
        """
        This needs to be updated to match configuration patterns used elsewhere in engine.  See #1904 for details.
        """
        self.center = definition['center']
        self.width = definition['width']
        self.strength = definition['strength']
        self.profile = definition['profile']
        self.emission_or_absorption = definition['emission_or_absorption']

        # In wavelength units
        self.wave_width = self.width * 1e3 / cs.c * self.center

        # This is to remind us of the internal units. An update could be to make this flexible, as it is in pysynphot.
        self.wunit = 'micron'
        self.vel_unit = 'km/s'
        self.funit = 'mJy'

    def line_waveset(self, wave, window_factor=5):
        """
        Create set of wavelengths to optimally sample the spectral line.

        Arguments
        ---------
        wave: np.ndarray
            Wavelength set of spectrum to which line will be added
        window_factor: int
            Configure the size of the window over which wavelength samples are checked and, if necessary, created.
            Window size is 2 * window_factor * width.

        Returns
        -------
        lwave: np.ndarray
            Array of wavelengths
        """
        nsamp = window_factor * 2 * 5  # Nyquist plus some oversampling
        wmin = np.max(self.center - self.wave_width * window_factor, 0.0)
        wmax = self.center + self.wave_width * window_factor

        # check input wavelength set over the line window to see if we need supplement with some more samples
        lsubs = (wave > wmin) & (wave < wmax)
        if len(wave[lsubs]) > nsamp:
            # if we're already sampled sufficiently over the line region, return wave unmodified
            lwave = wave
        else:
            # otherwise create a finely sampled set over the line, merge with the input, and return result
            lwaveset = np.linspace(wmin, wmax, int(nsamp))
            lwave = merge_wavelengths(wave, lwaveset)

        return lwave

    def flux(self, wave):
        """
        Calculate emission line flux over wavelength set, wave.

        Arguments
        ---------
        wave: np.ndarray
            Wavelengths over which to calculate line emission

        Returns
        -------
        flux: np.ndarray
            Line flux at wave
        """
        sigma = self.wave_width / 2.3548

        # normalized to a peak flux density of 1 mJy
        flux = np.exp(-(wave - self.center) ** 2 / (2. * sigma ** 2))

        # input line strength units in erg/cm^2/s
        int_flux = -ig.simps(flux * 1e-26, cs.c / (wave * 1e-6))

        flux = flux * self.strength / int_flux
        return flux

    def optical_depth(self, wave):
        """
        Calculate absorption line optical depth over wavelength set, wave.

        Arguments
        ---------
        wave: np.ndarray
            Wavelengths over which to calculate line optical depth

        Returns
        -------
        tau: np.ndarray
            Line optical depth at wave
        """
        sigma = self.wave_width / 2.3548
        tau = np.exp(-(wave - self.center) ** 2 / (2. * sigma ** 2))
        tau *= self.strength
        return tau


class ConvolvedSceneCube(object):

    """
    The SceneCube contains the source flux distribution as seen through the optics
    of a telescope/instrument combination. It takes the sampled ModelSceneCube and
    convolves it with the appropriate PSF from PSFLibrary.

    This is a central class of the ETC. The cube has two spatial and one wavelength dimension.

    Parameters
    ----------
    Sources : list
        A list of Source instances containing the physical parameters
        of the sources within the scene.
    instrument : Instrument class
        An instance of the instrument class
    background : Background, optional
        A background spectrum.
    PSFLibrary : PSFLibrary
        A library of the PSFs to use.


    Attributes
    ----------
    PSFLibrary :
    instrument :
    Grid :
    flux_cube :
    flux_plus_bg :

    """

    def __init__(self, scene, instrument, background=None, psf_library=None, webapp=False):
        self.warnings = {}
        self.scene = scene
        self.psf_library = psf_library
        self.aper_width = instrument.get_aperture_pars()['disp']
        self.aper_height = instrument.get_aperture_pars()['xdisp']
        self.multishutter = instrument.get_aperture_pars()['multishutter']
        nslice_str = instrument.get_aperture_pars()['nslice']

        if nslice_str is not None:
            self.nslice = int(nslice_str)
        else:
            self.nslice = 1

        self.instrument = instrument
        self.background = background

        self.fov_size = self.get_fov_size()

        # Figure out what the relevant wavelength range is, given the instrument mode
        wrange = self.current_instrument.get_wave_range()

        self.source_spectra = []

        # run through the sources and check their wavelength extents. warn if they fall short of the
        # current instrument configuration's range.
        mins = []
        maxes = []
        key = None
        for i, src in enumerate(scene.sources):
            spectrum = AstroSpectrum(src, webapp=webapp)
            self.warnings.update(spectrum.warnings)
            smin = spectrum.wave.min()
            smax = spectrum.wave.max()
            if smin > wrange['wmin']:
                if smin > wrange['wmax']:
                    key = "spectrum_missing_red"
                    msg = warning_messages[key] % (smin, smax, wrange['wmax'])
                    self.warnings["%s_%s" % (key, i)] = msg
                else:
                    key = "wavelength_truncated_blue"
                    msg = warning_messages[key] % (smin, wrange['wmin'])
                    self.warnings["%s_%s" % (key, i)] = msg
            if smax < wrange['wmax']:
                if smax < wrange['wmin']:
                    key = "spectrum_missing_blue"
                    msg = warning_messages[key] % (smin, smax, wrange['wmin'])
                    self.warnings["%s_%s" % (key, i)] = msg
                else:
                    key = "wavelength_truncated_red"
                    msg = warning_messages[key] % (smax, wrange['wmax'])
                    self.warnings["%s_%s" % (key, i)] = msg

            mins.append(smin)
            maxes.append(smax)

        wmin = max([np.array(mins).min(), wrange['wmin']])
        wmax = min([np.array(maxes).max(), wrange['wmax']])

        # make sure we have something within range and error out otherwise
        if wmax < wrange['wmin'] or wmin > wrange['wmax']:
            msg = "No wavelength overlap between source_spectra [%.2f, %.2f] and instrument [%.2f, %.2f]." % (
                wmin,
                wmax,
                wrange['wmin'],
                wrange['wmax']
            )
            raise RangeError(value=msg)

        # warn if partial overlap between combined wavelength range of all sources and the instrument's wrange
        if wmin != wrange['wmin'] or wmax != wrange['wmax']:
            key = "scene_range_truncated"
            self.warnings[key] = warning_messages[key] % (wmin, wmax, wrange['wmin'], wrange['wmax'])

        """
        Trim spectrum and do the spectral convolution here on a per-spectrum basis.  Most efficient to do it here
        before the wavelength sets are merged.  Also easier and much more efficient than convolving
        an axis of a 3D cube.
        """
        for src in scene.sources:
            spectrum = AstroSpectrum(src, webapp=webapp)
            # we trim here as an optimization so that we only convolve the section we need of a possibly very large spectrum
            spectrum.trim(wrange['wmin'], wrange['wmax'])
            spectrum = instrument.spectrometer_convolve(spectrum)
            self.source_spectra.append(spectrum)

        """
        different spectra will have different sets of wavelengths. the obvious future
        case will be user-supplied spectra, but this is also true for analytic spectra
        that have different emission/absorption lines. go through each of the spectra,
        merge all of the wavelengths sets into one, and then resample each spectrum
        onto the combined wavelength set.
        """
        self.wave = self.source_spectra[0].wave
        for s in self.source_spectra:
            self.wave = merge_wavelengths(self.wave, s.wave)

        projection_type = instrument.projection_type

        """
        For the spectral projections, we could use the pixel sampling. However, this
        may oversample the cube for input spectra with no narrow features. So we first check
        whether the pixel sampling will give us a speed advantage. Otherwise, do not resample to
        an unnecessarily fine wavelength grid.
        """
        if projection_type in ('spec', 'slitless', 'multiorder'):
            wave_pix = instrument.get_wave_pix()
            if wave_pix.size < self.wave.size:
                self.wave = wave_pix[np.where(np.logical_and(wave_pix >= wrange['wmin'],
                                                             wave_pix <= wrange['wmax']))]

        """
        There is no inherently optimal sampling for imaging modes, but we resample here to
        a reasonable number of wavelength bins if necessary. This helps keep the cube rendering reasonable
        for large input spectra. Note that the spectrum resampling now uses the flux conserving method
        of pysynphot.
        """
        if projection_type == 'image':
            """
            In practice a value of 200 samples within an imaging configuration's wavelength range
            (i.e. filter bandpass) should be more than enough. Note that because we use pysynphot
            to resample, the flux of even narrow lines is conserved.
            """
            nw_maximal = 200
            if self.wave.size > nw_maximal:
                self.wave = np.linspace(wrange['wmin'], wrange['wmax'], nw_maximal)

        self.nw = self.wave.size
        self.total_flux = np.zeros(self.nw)

        for spectrum in self.source_spectra:
            spectrum.resample(self.wave)
            self.total_flux += spectrum.flux

        # also need to resample the background spectrum
        if self.background is not None:
            self.background.resample(self.wave)

        self.grid, self.aperture_list, self.flux_cube_list, self.flux_plus_bg_list = \
            self.create_flux_cube(background=self.background)
        self.dist = self.grid.dist()

    def get_fov_size(self, pixbuffer=20):
        # The scene size is the minimum size containing all sources, but at least as large as the PSF.

        instrument_name = self.instrument.get_name()
        aperture_name = self.instrument.get_aperture()
        psf_shape = self.psf_library.get_shape(instrument_name, aperture_name)
        psf_pix_scl = self.psf_library.get_pix_scale(instrument_name, aperture_name)
        psf_size = psf_shape[0] * psf_pix_scl

        if self.instrument.projection_type == 'multiorder':
            # for the multiorder case (i.e. SOSS), we need to set the FOV size based on the subarray used
            subarray = self.instrument.detector['subarray']
            aperture = self.instrument.instrument['aperture']
            subarray_config = self.instrument.subarray_config[subarray]
            nx = subarray_config['nx']
            ny = subarray_config['ny']
            pix_scale = self.instrument.aperture_config[aperture]['pix']
            disp_axis = self.instrument.dispersion_axis()
            if disp_axis == "x":
                fov_pix = ny
            else:
                fov_pix = nx
            fov_size = fov_pix * pix_scale
        elif self.instrument.dynamic_scene:
            scene_size = self.scene.get_size()
            # a scene size and maximum scene size must be defined for each instrument/mode.
            # if they're not defined at all, it's a data problem.
            try:
                inst_fov = self.instrument.scene_size
                max_fov = self.instrument.max_scene_size
            except AttributeError as e:
                message = "Instrument configuration must specify default and maximum scene sizes. (%s)" % e
                raise DataError(value=message)
            # the configured scene size must be smaller than the maximum. since they can come from
            # either the data files or input configuration, this is an input error.
            if inst_fov > max_fov:
                message = "Specified scene size, %f arcsec, larger than maximum allowed size of %f arcsec" % (inst_fov, max_fov)
                raise EngineInputError(value=message)

            # find the largest of the size of the defined scene, the configured instrument FOV, and the PSF image size.
            # if this is larger than the configured maximum FOV size, set it to the maximum value. should add logic
            # here or in Report to emit a warning when this happens.
            fov_size = np.max([scene_size + pixbuffer * psf_pix_scl, inst_fov, psf_size])
            if fov_size > max_fov:
                key = "max_scene_size_reached"
                self.warnings[key] = warning_messages[key] % (fov_size, max_fov)
                fov_size = max_fov
        else:
            fov_size = psf_size

        # check to make sure at least one source is within the field of view. warn otherwise...
        if fov_size < self.scene.get_min_size():
            key = "scene_fov_too_small"
            self.warnings[key] = warning_messages[key] % fov_size

        return fov_size

    def create_flux_cube(self, background=None):
        """
        Generate the list of convolved flux cubes that will go into the ETC calculation.
        The spectral convolution is already done and the spatial convolution is performed here
        using self.PSFLibrary (nominally as generated by webbPSF). The flux cube handles position-dependent
        PSFs by assigning PSF profiles to individual sources, convolving intermediate cubes of sources with
        common PSFs and finally co-adding cubes to create a final master flux cube. Typical ETC calculations, which do
        not have position-dependent PSFs are not affected by this functionality.

        Parameters
        ----------
        background: background.Background instance

        Returns
        -------
        <tuple>:
            spatial grid used to create cube(s) (coords.Grid instance)
            list of apertures (list)
            list of flux cubes (list; one per aperture)
            list of flux cuves including background (list; one per aperture)

        """
        if self.psf_library is None:
            raise NotImplementedError('The use of a simple PSF in isolation is deprecated. Please provide PSFLibrary.')

        instrument_name = self.instrument.get_name()
        aperture_name = self.instrument.get_aperture()
        psf_pixsize = self.psf_library.get_pix_scale(instrument_name, aperture_name)
        psf_upsamp = self.psf_library.get_upsamp(instrument_name, aperture_name)

        detector_npix = int(np.round(self.fov_size / psf_pixsize / psf_upsamp))
        if detector_npix % 2 == 0:
            detector_npix += 1
        detector_shape = (detector_npix, detector_npix)

        flux_cube_list = [
            np.zeros(
                (detector_shape[0],
                 detector_shape[1],
                 self.nw)) for ir in range(self.nslice)]

        flux_plus_bg_list = [
            np.zeros(
                (detector_shape[0],
                 detector_shape[1],
                 self.nw)) for ir in range(self.nslice)]

        if background is not None:
            self.bg = background.mjy_pix
        else:
            self.bg = self.wave * 0.0

        scene_grid = Grid(psf_pixsize, psf_pixsize, detector_shape[0] * psf_upsamp, detector_shape[1] * psf_upsamp)

        psf_associations = self.psf_library.associate_offset_to_source(self.scene.sources, instrument_name, aperture_name)
        unique_offsets = list(set(psf_associations))

        current_scenes = []
        for unique_offset in unique_offsets:
            offset_indices = [i for (i, v) in enumerate(psf_associations) if v == unique_offset]
            current_scene = ModelSceneCube([self.source_spectra[i] for i in offset_indices], scene_grid)
            current_scenes.append(current_scene)

        for iw in np.arange(self.nw):
            for current_scene, unique_offset, i in zip(current_scenes, unique_offsets, range(len(unique_offsets))):
                if i == 0:
                    psf = AdvancedPSF(
                        self.wave[iw],
                        self.psf_library,
                        instrument_name,
                        aperture_name,
                        psf_source_offset=unique_offset,
                        model_scene=current_scene,
                        aper_width=self.aper_width,
                        aper_height=self.aper_height,
                        multishutter=self.multishutter,
                        nslice=self.nslice,
                        bg_w=self.bg[iw]
                    )
                else:
                    # if there are sources with different PSFs, calculate their intensities and add them,
                    # but do not add more background.
                    psf.add_intensity(
                        AdvancedPSF(
                            self.wave[iw],
                            self.psf_library,
                            instrument_name,
                            aperture_name,
                            psf_source_offset=unique_offset,
                            model_scene=current_scene,
                            aper_width=self.aper_width,
                            aper_height=self.aper_height,
                            multishutter=self.multishutter,
                            nslice=self.nslice,
                            bg_w=0.
                        )
                    )

            for islice in range(self.nslice):
                flux_cube_list[islice][:, :, iw] = psf.slice_int_list[islice]
                flux_plus_bg_list[islice][:, :, iw] = psf.slice_int_plus_bg_list[islice]

        return psf.grid, psf.aperture_list, flux_cube_list, flux_plus_bg_list

    def spectral_model_transform(self):
        """
        Create engine API format dict section containing properties of the wavelength coordinates
        used in the construction of a ConvolvedSceneCube.

        Returns
        -------
        t: dict (engine API compliant keys)
        """
        t = {}
        t['wave_refpix'] = 0
        t['wave_refval'] = self.wave[0]
        t['wave_max'] = self.wave.max()
        t['wave_min'] = self.wave.min()
        t['wave_size'] = self.wave.size
        # this is a bit of a hack since the wavelength sampling is NOT constant at this stage.
        # this is simply the mean step and should probably be removed altogether. if one wishes to
        # plot anything from this stage, they should instead use self.wave which provides the true
        # mapping of index -> wavelength
        if len(self.wave > 1):
            t['wave_step'] = (self.wave.max() - self.wave.min()) / self.wave.size
        else:
            t['wave_step'] = 0.0
        return t

    def cube_wcs_info(self):
        """
        Create WCS headers and FITS binary table that describe the cube's coordinate system.

        Returns
        -------
        tbhdu: astropy.io.fits.BinTableHDU instance
            The wavelength sampling is irregularly spaced so we define a binary FITS
            table that contains the array of wavelengths, self.wave.
        header: dict
            Contains the WCS keys that define the coordinate transformation for all axes
        """
        header = self.grid.wcs_info()
        header['ctype3'] = 'WAVE-TAB'
        header['cname3'] = 'Wavelength'
        header['cunit3'] = 'um'
        header['PS3_0'] = 'WCS-TAB'
        header['PS3_1'] = 'WAVELENGTH'
        header['PS3_2'] = 'WAVE-INDEX'
        tbhdu = fits.BinTableHDU.from_columns([
            fits.Column(name='WAVELENGTH',
                        unit='um',
                        format="1D",
                        array=self.wave),
            fits.Column(name='WAVE-INDEX',
                        format="1J",
                        array=np.arange(self.wave.size))
        ])
        tbhdu.name = 'WCS-TAB'
        return tbhdu, header

    def export_to_fits(self, fitsfile='ModelDetectorCube'):
        header = self.wcs_info()
        slice_index = 0
        for flux_plus_bg in self.flux_plus_bg_list:
            fitsfile_slice = fitsfile + str(slice_index).strip() + '.fits'
            fits.writeto(fitsfile_slice, np.rollaxis(flux_plus_bg, 2), header, clobber=True)
            fits.append(fitsfile_slice, self.wave)


class AdvancedPSF(object):

    """
    Convolve scene with PSF from a PSF library (e.g., one calculated using WebbPSF).

    Parameters
    ----------
    wave: float
        Wavelength (in microns) at which the PSF is to be calculated
    psf_library: psf_library.PSFLibrary instance
        Library of PSFs to use for convolution kernels
    instrument: string
        Name of instrument being used
    mode: string
        Name of observing mode being used
    model_scene: astro_spectrum.ModelSceneCube instance
        Model scene as sampled into an image cube
    aper_width: float
        Width of the aperture used to slice the PSF
    aper_height: float
        Height of the aperture used to slice the PSF
    multishutter: list of (float, float)
        List of tuples to construct slit apertures as a set of identical rectangles.
        Each tuple in the list is the X and Y offset of a rectangle.
    nslice: int
        Number of slices
    bg_w: float
        Background rate at wavelength, wave.
    """

    def __init__(self, wave, psf_library, instrument, mode, model_scene=None, psf_source_offset=(0, 0),
                 aper_width=None, aper_height=None, multishutter=[(0.0, 0.0)], nslice=None, bg_w=0):

        # The advanced PSF gets its intensity map from a PSF library
        profile = psf_library.get_psf(wave, instrument, mode, source_offset=psf_source_offset)
        psf_upsamp = profile['upsamp']
        psf_pixscl = profile['pix_scl']

        if profile['int'].shape[0] != profile['int'].shape[1]:
            raise ValueError("The PSF must have a square grid shape nx=ny")
#        psf_npix = profile['int'].shape[0]
        scene_npix = model_scene.nx

        kernel = profile['int']
#        npix = psf_npix / psf_upsamp
        npix = scene_npix / psf_upsamp
        new_shape = (npix, npix)

        # If there is a scene, convolve with it
        if model_scene is not None:
            # close enough?
            if (np.abs(model_scene.xsamp / psf_pixscl - 1) > 1e-10):
                raise ValueError("scene sampling must be the same as PSF sampling")
            windex = self._find_nearest_index(wave, model_scene.wave)
            """
            See https://github.com/STScI-SSB/scamp/issues/71 for discussion and benchmarks about
            these different FFT methods. astropy.convolution plus pyFFTW yields a significant
            speed-up, of order 30% or more.
            """
            # The original method
            # import stsci.convolve as stsci
            # self.intensity = stsci.convolve2d(model_scene.int[:, :, windex], kernel, fft=1, mode='wrap')

            # Basic astropy.convolution with FFTs with parameters set to produce numbers that
            # match the original method
            # self.intensity = convolve_fft(model_scene.int[:, :, windex], kernel[:-1, :-1], normalize_kernel=False)

            # Optimal solution using astropy.convolution plus pyFFTW
            self.intensity = convolve_fft(model_scene.int[:, :, windex], kernel, normalize_kernel=False,
                                          boundary='fill',
                                          fftn=pyfftw.interfaces.numpy_fft.fftn,
                                          ifftn=pyfftw.interfaces.numpy_fft.ifftn)
        else:
            self.intensity = profile['int']

        #
        # Add the background
        self.intensity_plus_bg = self.intensity + bg_w / psf_upsamp ** 2

        """
        We can now operate with any number of physical spectral apertures (slices) of the FOV. A single slit
        mode simply has nslice=1. An imaging mode is also a slice, but with infinite aperture.
        A multishutter instrument can create a slice aperture mask consisting of a discrete number of mutually
        offset rectangles. In principle, one could create an IFU with each slice consisting of discrete shutters.
        This could be used to simulate different IFU designs, such as lenslet or micro-mirror arrays.
        """
        if aper_width is not None and aper_height is not None:
            offsets = [(i - (nslice - 1) / 2.) * aper_width for i in np.arange(nslice)]
            self.slice_int_list = []
            self.slice_int_plus_bg_list = []
            self.slice_mask_list = []
            self.aperture_list = []
            self.grid = Grid(psf_pixscl * psf_upsamp, psf_pixscl * psf_upsamp, npix, npix)
            fine_grid = Grid(psf_pixscl, psf_pixscl, scene_npix, scene_npix)
            for offset in offsets:
                slice_mask_fine = np.zeros(self.intensity.shape)
                # Is this a multishutter instrument?
                if multishutter:
                    for shutter in multishutter:
                        new_mask = fine_grid.rectangular_mask(
                            width=aper_width,
                            height=aper_height,
                            xoff=offset + shutter[0],
                            yoff=shutter[1]
                        )
                        slice_mask_fine = np.maximum(slice_mask_fine, new_mask)
                else:
                    slice_mask_fine = fine_grid.rectangular_mask(
                        width=aper_width,
                        height=aper_height,
                        xoff=offset,
                        yoff=0.0
                    )

                slice_int, slice_int_plus_bg, slice_mask = self._apply_slit_mask(slice_mask_fine, new_shape)
                # for this purpose a set of multiple shutters is treated as a single aperture and uses
                # the properties of the central shutter.
                aperture = {'width': aper_width, 'height': aper_height, 'offset': (0., offset)}
                self.slice_int_list.append(slice_int)
                self.slice_int_plus_bg_list.append(slice_int_plus_bg)
                self.slice_mask_list.append(slice_mask)
                self.aperture_list.append(aperture)

        else:
            slice_mask_fine = np.ones((scene_npix, scene_npix))
            slice_int, slice_int_plus_bg, slice_mask = self._apply_slit_mask(slice_mask_fine, new_shape)
            self.grid = Grid(psf_pixscl * psf_upsamp, psf_pixscl * psf_upsamp, npix, npix)
            self.slice_int_list = [slice_int]
            self.slice_int_plus_bg_list = [slice_int_plus_bg]
            self.slice_mask_list = [slice_mask]
            self.aperture_list = [self.grid.get_aperture()]

    def add_intensity(self, psf):
        """
        Add another compatible PSF intensity.

        Parameters
        ----------
        psf : AdvancedPSF instance
            A previously computed AdvancedPSF to add.

        """

        self.slice_int_list = [slice_int + add_slice_int for slice_int, add_slice_int in
                               zip(self.slice_int_list, psf.slice_int_list)]
        self.slice_int_plus_bg_list = [slice_int + add_slice_int for slice_int, add_slice_int in
                                       zip(self.slice_int_plus_bg_list, psf.slice_int_plus_bg_list)]

    def _find_nearest_index(self, number, vector):
        """
        Find index of 'vector' closest to 'number'

        Parameters
        ----------
        number : float-like
            Value to compare for finding nearest index
        vector : ndarray
            Vector to search for nearest index

        Returns
        -------
        float
        """
        index = np.abs(vector - number).argmin()
        return index

    def _apply_slit_mask(self, slit_mask, new_shape):
        """
        Applies a slit mask, and also rebins to a new shape.

        Parameters
        ----------
        slit_mask : 2D ndarray
            Slit mask as generated by self._create_slit_mask()
        new_shape : list-like
            New shape for array after binning

        Returns: list-like of 2D ndarrays
            Intensity image, intensity plus background image, aperture mask
        """
        intensity = self._rebin(self.intensity * slit_mask, new_shape)
        intensity_plus_bg = self._rebin(self.intensity_plus_bg * slit_mask, new_shape)
        aperture_mask = self._rebin(slit_mask, new_shape)
        return intensity, intensity_plus_bg, aperture_mask

    def _rebin(self, a, shape):
        """
        Re-bin a 2D array.

        Parameters
        ----------
        a : ndarray
            Array to be re-binned
        shape : list-like
            New shape after re-binning

        Returns
        -------
        ndarray
        """
        sh = int(shape[0]), int(a.shape[0] // shape[0]), int(shape[1]), int(a.shape[1] // shape[1])
        new = a.reshape(sh).sum(-1).sum(1)
        return new

    def _rebin_1d_mean(self, a, shape):
        """
        Re-bin a 1D array using averaging.

        Parameters
        ----------
        a : 1D ndarray
            Array to be re-binned
        shape : int
            New length after re-binning

        Returns
        -------
        ndarray
        """
        sh = shape, a.shape[1] // shape
        new = a.reshape(sh).mean(-1)
        return new
