# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import copy
import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as sci_int
from scipy.ndimage.interpolation import shift
from astropy.convolution import convolve_fft
import pyfftw

from . import observation
from . import astro_spectrum as astro
from . import background as bg
from . import coords
from .config import DefaultConfig
from .report import Report
from .scene import Scene
from .calc_utils import build_empty_scene
from .custom_exceptions import EngineInputError, EngineOutputError, RangeError, DataError
from .instrument_factory import InstrumentFactory
from .strategy import StrategyFactory

from six.moves import zip
# pyfftw.interfaces.cache.enable()


class CalculationConfig(DefaultConfig):

    """
    Encapsulate calculation configuration parameters (e.g. effects to include, noise sources to consider)
    """
    pass


class DetectorSignal(astro.ConvolvedSceneCube):

    """
    This class contains functionality for calculating the integrated electron rate of an
    astronomical source at the relevant instrument detector plane.

    Parameters
    ----------
    observation: observation.Observation instance
        Contains information required to configure a calculation
    calc_config: CalculationConfig instance
        Contains boolean flags that control which noise components are included in the calculation
    webapp: bool
        Toggle strict engine API checking
    order: int or None
        For multi-order case (i.e. SOSS), designate which order to use

    Attributes
    ----------
    """

    def __init__(self, observation, calc_config=CalculationConfig(), webapp=False, order=None):
        # Get calculation configuration
        self.calculation_config = calc_config

        # Link to the passed observation
        self.observation = observation

        # Load the instrument we're using
        self.current_instrument = observation.instrument
        # and configure it for the order we wish to use, if applicable
        self.current_instrument.order = order

        # how are we projecting the signal onto the detector plane?
        self.projection_type = self.current_instrument.projection_type

        # If we're in a dispersed mode, we need to know which axis the signal is dispersed along
        self.dispersion_axis = self.current_instrument.dispersion_axis()

        # Get the detector parameters (read noise, etc.)
        self.det_pars = self.current_instrument.get_detector_pars()

        # Initialize detector mask
        self.det_mask = 1.0

        # Get the background
        if self.calculation_config.effects['background']:
            self.background = bg.Background(self.observation, webapp=webapp)
        else:
            self.background = None

        # Then initialize the flux and wavelength grid
        astro.ConvolvedSceneCube.__init__(
            self,
            self.observation.scene,
            self.current_instrument,
            background=self.background,
            psf_library=self.current_instrument.psf_library,
            webapp=webapp
        )

        self.warnings.update(self.background.warnings)

        # We have to propagate the background through the system transmission
        # to get the background in e-/s/pixel/micron. The background rate is a 1D spectrum.
        self.bg_fp_rate = self.get_bg_fp_rate()

        # Initialize slice lists
        self.rate_list = []
        self.rate_plus_bg_list = []
        self.saturation_list = []
        self.pixgrid_list = []

        # Loop over all slices and calculate the photon and electron rates through the
        # observatory for each one. Note that many modes (imaging, etc.) will have just
        # a single slice.
        for flux_cube, flux_plus_bg in zip(self.flux_cube_list, self.flux_plus_bg_list):
            # Rates for the slice without the background
            slice_rate = self.all_rates(flux_cube, add_extended_background=False)

            # Rates for the slice with the background added
            slice_rate_plus_bg = self.all_rates(flux_plus_bg, add_extended_background=True)

            # Saturation map for the slice
            slice_saturation = self.get_saturation_mask(rate=slice_rate_plus_bg['fp_pix'])

            # The grid in the slice
            slice_pixgrid = self.get_pix_grid(slice_rate)

            # Append all slices to the master lists
            self.rate_list.append(slice_rate)
            self.rate_plus_bg_list.append(slice_rate_plus_bg)
            self.saturation_list.append(slice_saturation)
            self.pixgrid_list.append(slice_pixgrid)

        # Get the mapping of wavelength to pixels on the detector plane. This is grabbed from the
        # first entry in self.rate_list and is currently defined to be the same for all slices.
        self.wave_pix = self.get_wave_pix()

        # This is also grabbed from the first slice as a diagnostic
        self.fp_rate = self.get_fp_rate()

        # Note that the 2D image due to background alone may have spatial structure due to instrumental effects.
        # Therefore it is calculated here.
        self.bg_pix_rate = self.get_bg_pix_rate()

        # Reassemble rates of multiple slices on the detector
        self.rate = self.on_detector(self.rate_list)
        self.rate_plus_bg = self.on_detector(self.rate_plus_bg_list)

        self.detector_pixels = self.current_instrument.get_detector_pixels(self.wave_pix)

        # Get the read noise correlation matrix and store it as an attribute.
        if self.det_pars['rn_correlation']:
            self.read_noise_correlation_matrix = self.current_instrument.get_readnoise_correlation_matrix(self.rate.shape)

    def spectral_detector_transform(self):
        """
        Create engine API format dict section containing properties of wavelength coordinates
        at the detector plane.

        Returns
        -------
        t: dict (engine API compliant keys)
        """
        t = {}
        t['wave_det_refpix'] = 0
        t['wave_det_max'] = self.wave_pix.max()
        t['wave_det_min'] = self.wave_pix.min()

        # there are currently three projection_type's which are basically detector plane types:
        #
        # 'spec' - where the detector plane is purely dispersion vs. spatial
        # 'slitless' - basically a special case of 'spec' with where dispersion and spatial are mixed
        # 'image' - where the detector plane is purely spatial vs. spatial (i.e. no disperser element)
        #
        # 'IFU' mode is of projection_type='spec' because the mapping from detector X pixels to
        # wavelength is the same for each slice.  this projection_type will work for 'MSA' mode as well
        # because we will only handle one aperture at a time.  'slitless' spectroscopy will mix
        # spatial and dispersion information onto the detector X axis.  however, the detector
        # plane is fundamentally spatial vs. wavelength in that case so it's handled the same as
        # projection_type='spec'. creating a spectrum for a specific target will be handled via the
        # extraction strategy.
        if self.projection_type in ('spec', 'slitless', 'multiorder'):
            t['wave_det_size'] = len(self.wave_pix)
            if len(self.wave_pix) > 1:
                # we don't yet have a way of handling non-linear coordinate transforms here. that said,
                # this is mostly right for most of our cases with nirspec prism being the notable exception.
                # this is also only used for plotting purposes while the true actual wave_pix mapping is used
                # internally for all calculations.
                t['wave_det_step'] = (self.wave_pix[-1] - self.wave_pix[0]) / t['wave_det_size']
            else:
                t['wave_det_step'] = 0.0
            t['wave_det_refval'] = self.wave_pix[0]
        elif self.projection_type == "image":
            t['wave_det_step'] = 0.0
            t['wave_det_refval'] = self.wave_pix[0]
            t['wave_det_size'] = 1
        else:
            message = "Unsupported projection_type: %s" % self.projection_type
            raise EngineOutputError(value=message)
        return t

    def wcs_info(self):
        """
        Get detector coordinate transform as a dict of WCS keyword/value pairs.

        Returns
        -------
        header: dict
            WCS header keys defining coordinate transform in the detector plane
        """
        if self.projection_type == 'image':
            # if we're in imaging mode, the detector sampling is the same as the model
            header = self.grid.wcs_info()
        elif self.projection_type in ('spec', 'slitless', 'multiorder'):
            # if we're in a dispersed mode, dispersion can be either along the X or Y axis. the image outputs in
            # the engine Report are rotated so that dispersion will always appear to be along the X axis with
            # wavelength increasing with increasing X (i. e. dispersion angle of 0).  currently, the only other
            # supported dispersion angle is 90 which is what we get when dispersion_axis == 'y'.
            t = self.grid.as_dict()
            t.update(self.spectral_detector_transform())
            header = {
                'ctype1': 'Wavelength',
                'crpix1': 1,
                'crval1': t['wave_det_min'] - 0.5 * t['wave_det_step'],
                'cdelt1': t['wave_det_step'],
                'cunit1': 'um',
                'cname1': 'Wavelength',
                'ctype2': 'Y offset',
                'crpix2': 1,
                'crval2': t['y_min'] - 0.5 * t['y_step'],
                'cdelt2': -t['y_step'],
                'cunit2': 'arcsec',
                'cname2': 'Detector Offset',
            }
            if self.dispersion_axis == 'y':
                header['ctype2'] = 'X offset'
                header['crval2'] = t['x_min'] - 0.5 * t['x_step'],
                header['cdelt2'] = t['x_step']
        else:
            message = "Unsupported projection_type: %s" % self.projection_type
            raise EngineOutputError(value=message)
        return header

    def get_wave_pix(self):
        """
        Return the mapping of wavelengths to pixels on the detector plane
        """
        return self.rate_list[0]['wave_pix']

    def get_fp_rate(self):
        """
        Return scene flux at the focal plane in e-/s/pixel/micron (excludes background)
        """
        return self.rate_list[0]['fp']

    def get_bg_fp_rate(self):
        """
        Calculate background in e-/s/pixel/micron at the focal plane
        """
        bg_fp_rate = self.focal_plane_rate(self.ote_rate(self.background.mjy_pix))
        return bg_fp_rate

    def get_bg_pix_rate(self):
        """
        Calculate the background on the detector in e-/s/pixel
        """
        bg_pix_rate = self.rate_plus_bg_list[0]['fp_pix'] - self.rate_list[0]['fp_pix']
        return bg_pix_rate

    def on_detector(self, rate_list):
        """
        This will take the list of (pixel) rates and use them create a single detector frame. A single
        image will only have one rate in the list, but the IFUs will have n_slices. There may be other examples,
        such as different spectral orders for NIRISS. It is not yet clear how many different flavors there are, so
        this step may get refactored if it gets too complicated. Observing modes that only have one set of rates
        (imaging and single-slit spectroscopy, for instance) will still go through this, but the operation is trivial.
        """
        aperture_sh = rate_list[0]['fp_pix'].shape
        n_apertures = len(rate_list)
        detector_shape = (aperture_sh[0] * n_apertures, aperture_sh[1])
        detector = np.zeros(detector_shape)

        i = 0
        for rate in rate_list:
            detector[i * aperture_sh[0]:(i + 1) * aperture_sh[0], :] = rate['fp_pix']
            i += 1

        return detector

    def get_pix_grid(self, rate):
        """
        Generate the coordinate grid of the detector plane
        """
        if self.projection_type == 'image':
            grid = self.grid
        elif self.projection_type in ('spec', 'slitless', 'multiorder'):
            nw = rate['wave_pix'].shape[0]
            if self.dispersion_axis == 'x':
                # for slitless calculations, the dispersion axis is longer than the spectrum being dispersed
                # because the whole field of view is being dispersed. 'excess' is the size of the FOV
                # and half will be to the left of the blue end of the spectrum and half to the right of the red end.
                # this is used to create the new spatial coordinate transform for the pixel image on the detector.
                excess = rate['fp_pix'].shape[1] - nw
                pix_grid = coords.IrregularGrid(
                    self.grid.col,
                    (np.arange(nw + excess) - (nw + excess) / 2.0) * self.grid.xsamp
                )
            else:
                excess = rate['fp_pix'].shape[0] - nw
                pix_grid = coords.IrregularGrid(
                    (np.arange(nw + excess) - (nw + excess) / 2.0) * self.grid.ysamp,
                    self.grid.row
                )
            return pix_grid
        else:
            raise EngineOutputError(value="Unsupported projection_type: %s" % self.projection_type)
        return grid

    def all_rates(self, flux, add_extended_background=False):
        """
        Calculate rates in e-/s/pixel/micron or e-/s/pixel given a flux cube in mJy

        Parameters
        ----------
        flux: ConvolvedSceneCube instance
            Convolved source flux cube with flux units in mJy
        add_extended_background: bool (default=False)
            Toggle for including extended background not contained within the flux cube

        Returns
        -------
        products: dict
            Dict of products produced by rate calculation.
                'wave_pix' - Mapping of wavelength to detector pixels
                'ote' - Source rate at the telescope aperture
                'fp' - Source rate at the focal plane in e-/s/pixel/micron
                'fp_pix' - Source rate per pixel
                'fp_pix_no_ipc' - Source rate per pixel excluding effects if inter-pixel capacitance
        """
        # The source rate at the telescope aperture
        ote_rate = self.ote_rate(flux)

        # The source rate at the focal plane in interacting photons/s/pixel/micron
        fp_rate = self.focal_plane_rate(ote_rate)

        # the fp_pix_variance is the variance of the per-pixel electron rate and includes the chromatic effects
        # of quantum yield.
        if self.projection_type == 'image':
            # The wavelength-integrated rate in e-/s/pixel, relevant for imagers
            fp_pix_rate, fp_pix_variance = self.image_rate(fp_rate)
            wave_pix = self.wave_eff(fp_rate)

        elif self.projection_type == 'spec':
            # The wavelength-integrated rate in e-/s/pixel, relevant for spectroscopy
            wave_pix, fp_pix_rate, fp_pix_variance = self.spec_rate(fp_rate)

        elif self.projection_type in ('slitless', 'multiorder'):
            # The wavelength-integrated rate in e-/s/pixel, relevant for slitless spectroscopy
            wave_pix, fp_pix_rate, fp_pix_variance = self.slitless_rate(
                fp_rate,
                add_extended_background=add_extended_background
            )

        else:
            raise EngineOutputError(value="Unsupported projection_type: %s" % self.projection_type)

        # Include IPC effects, if available and requested
        if self.det_pars['ipc'] and self.calculation_config.effects['ipc']:
            kernel = self.current_instrument.get_ipc_kernel()
            fp_pix_rate_ipc = self.ipc_convolve(fp_pix_rate, kernel)
        else:
            fp_pix_rate_ipc = fp_pix_rate

        # fp_pix is the final product. Since there is no reason to
        # carry around the ipc label everywhere, we rename it here.
        products = {
            'wave_pix': wave_pix,
            'ote': ote_rate,
            'fp': fp_rate,
            'fp_pix': fp_pix_rate_ipc,
            'fp_pix_no_ipc': fp_pix_rate,  # this is for calculating saturation
            'fp_pix_variance': fp_pix_variance  # this is for calculating the detector noise
        }

        return products

    def ote_rate(self, flux):
        """
        Calculate source rate in e-/s/pixel/micron at the telescope entrance aperture given
        a flux cube in mJy/pixel.
        """
        # spectrum in mJy/pixel, wave in micron, f_lambda in photons/cm^2/s/micron
        f_lambda = 1.5091905 * (flux / self.wave)
        ote_int = self.current_instrument.telescope.get_ote_eff(self.wave)
        coll_area = self.current_instrument.telescope.coll_area
        a_lambda = coll_area * ote_int
        # e-/s/pixel/micron
        ote_rate = f_lambda * a_lambda
        return ote_rate

    def focal_plane_rate(self, rate):
        """
        Takes the output from self.ote_rate() and multiplies it by the components
        of efficiency within the system and returns the source rate at the focal plane in
        e-/s/pixel/micron.
        """
        filter_eff = self.current_instrument.get_filter_eff(self.wave)
        disperser_eff = self.current_instrument.get_disperser_eff(self.wave)
        internal_eff = self.current_instrument.get_internal_eff(self.wave)
        qe = self.current_instrument.get_detector_qe(self.wave)

        fp_rate = rate * filter_eff * disperser_eff * internal_eff * qe
        return fp_rate

    def spec_rate(self, rate):
        '''
        For slitted spectrographs, calculate the detector signal by integrating
        along the dispersion direction of the cube (which is masked by a, by assumption,
        narrow slit). For slitless systems or slits wider than the PSF, the slitless_rate
        method should be used to preserve spatial information within the slit.

        Parameters
        ---------
        rate: numpy.ndarray
            Rate of photons interacting with detector as a function of model wavelength set

        Returns
        -------
        products: 3-element tuple of numpy.ndarrays
            first element - map of pixel to wavelength
            second element - electron rate per pixel
            third element - variance of electron rate per pixel
        '''
        dispersion = self.current_instrument.get_dispersion(self.wave)
        wave_pix = self.current_instrument.get_wave_pix()
        wave_pix_trunc = wave_pix[np.where(np.logical_and(wave_pix >= self.wave.min(),
                                                          wave_pix <= self.wave.max()))]

        # Check that the source spectrum is actually inside the instrumental wavelength
        # coverage.
        if len(wave_pix_trunc) == 0:
            raise RangeError(value='wave and wave_pix do not overlap')

        # Check the dispersion axis to determine which axis to sum and interpolate over
        if self.dispersion_axis == 'x':
            axis = 1
        else:
            axis = 0

        # We can simply sum over the dispersion direction. This is where we lose the spatial information within the aperture.
        spec_rate = np.sum(rate, axis=axis)

        # And then scale to the dispersion function (pixel/micron) to transform
        # from e-/s/micron to e-/s/pixel.
        spec_rate_pix = spec_rate * dispersion

        # but we are still sampled on the internal grid, so we have to interpolate to the pixel grid.
        # use kind='slinear' since it's ~2x more memory efficient than 'linear'. 'slinear' uses different code path to
        # calculate the slopes.
        int_spec_rate = sci_int.interp1d(self.wave, spec_rate_pix, axis=axis, kind='slinear', assume_sorted=True, copy=False)
        spec_rate_pix_sampled = int_spec_rate(wave_pix_trunc)

        # Handle a detector gap here by constructing a mask. If the current_instrument implements it,
        # it'll be a real mask array.  Otherwise it will simply be 1.0.
        self.det_mask = self.current_instrument.create_gap_mask(wave_pix_trunc)

        # this is the interacting photon rate in the detector with mask applied.
        spec_rate_pix_sampled *= self.det_mask

        # Add effects of non-unity quantum yields. For the spec projection, we assume that the quantum yield does not
        # change over a spectral element. Then we can just multiply the products by the relevant factors.
        q_yield, fano_factor = self.current_instrument.get_quantum_yield(wave_pix_trunc)

        # convert the photon rate to electron rate by multiplying by the quantum yield which is a function of wavelength
        spec_electron_rate_pix = spec_rate_pix_sampled * q_yield

        # the variance in the electron rate, Ve, is also scaled by the quantum yield plus a fano factor which is
        # analytic in the simple 1 or 2 electron case: Ve = (qy + fano) * Re.  since Re is the photon rate
        # scaled by the quantum yield, Re = qy * Rp, we get: Ve = qy * (qy + fano) * Rp
        spec_electron_variance_pix = spec_rate_pix_sampled * q_yield * (q_yield + fano_factor)
        products = wave_pix_trunc, spec_electron_rate_pix, spec_electron_variance_pix

        return products

    def image_rate(self, rate):
        '''
        Calculate the electron rate for imaging modes by integrating along
        the wavelength direction of the cube.

        Parameters
        ---------
        rate: numpy.ndarray
            Rate of photons interacting with detector as a function of model wavelength set

        Returns
        -------
        products: 2-element tuple of numpy.ndarrays
            first element - electron rate per pixel
            second element - variance of electron rate per pixel
        '''
        q_yield, fano_factor = self.current_instrument.get_quantum_yield(self.wave)

        # convert the photon rate to electron rate by multiplying by the quantum yield which is a function of wavelength
        electron_rate_pix = integrate.simps(rate * q_yield, self.wave)

        # the variance in the electron rate, Ve, is also scaled by the quantum yield plus a fano factor which is
        # analytic in the simple 1 or 2 electron case: Ve = (qy + fano) * Re.  since Re is the photon rate
        # scaled by the quantum yield, Re = qy * Rp, we get: Ve = qy * (qy + fano) * Rp
        electron_variance_pix = integrate.simps(rate * q_yield * (q_yield + fano_factor), self.wave)

        products = electron_rate_pix, electron_variance_pix

        return products

    def slitless_rate(self, rate, add_extended_background=True):
        '''
        Calculate the detector rates for slitless modes. Here we retain all spatial information and build
        up the detector plane by shifting and coadding the frames from the convolved flux cube. Also need to handle
        and add background that comes from outside the flux cube, but needs to be accounted for.

        Parameters
        ----------
        rate: 3D numpy.ndarray
            Cube containing the flux rate at the focal plane
        add_extended_background: bool (default: True)
            Toggle for including extended background not contained within the flux cube

        Returns
        -------
        products: 2 entry tuple
            wave_pix: 1D numpy.ndarray containing wavelength to pixel mapping on the detector plane
            spec_rate: 2D numpy.ndarray of detector count rates
        '''
        wave_pix = self.current_instrument.get_wave_pix()
        wave_subs = np.where(
            np.logical_and(
                wave_pix >= self.wave.min(),
                wave_pix <= self.wave.max()
            )
        )
        wave_pix_trunc = wave_pix[wave_subs]

        if len(wave_pix_trunc) == 0:
            raise RangeError(value='wave and wave_pix do not overlap')

        dispersion = self.current_instrument.get_dispersion(wave_pix_trunc)
        trace = self.current_instrument.get_trace(wave_pix_trunc)

        q_yield, fano_factor = self.current_instrument.get_quantum_yield(wave_pix_trunc)

        # use kind='slinear' since it's ~2x more memory efficient than 'linear'. 'slinear' uses different code path to
        # calculate the slopes.
        int_rate_pix = sci_int.interp1d(self.wave, rate, kind='slinear', axis=2, assume_sorted=True, copy=False)
        rate_pix = int_rate_pix(wave_pix_trunc)

        # convert the photon rate to electron rate by multiplying by the quantum yield which is a function of wavelength
        electron_rate_pix = rate_pix * q_yield

        # the variance in the electron rate, Ve, is also scaled by the quantum yield plus a fano factor which is
        # analytic in the simple 1 or 2 electron case: Ve = (qy + fano) * Re.  since Re is the photon rate
        # scaled by the quantum yield, Re = qy * Rp, we get: Ve = qy * (qy + fano) * Rp
        electron_variance_pix = rate_pix * q_yield * (q_yield + fano_factor)

        # interpolate the background onto the pixel spacing
        int_bg_fp_rate = sci_int.interp1d(self.wave, self.bg_fp_rate, kind='slinear', assume_sorted=True, copy=False)
        bg_fp_rate_pix = int_bg_fp_rate(wave_pix_trunc)

        # calculate electron rate and variance due to background
        bg_electron_rate = bg_fp_rate_pix * q_yield
        bg_electron_variance = bg_fp_rate_pix * q_yield * (q_yield + fano_factor)

        # dispersion_axis tells us whether we need to sum the planes of the cube horizontally
        # or vertically on the detector plane.
        if self.dispersion_axis == 'x':
            spec_shape = (rate_pix.shape[0], rate_pix.shape[2] + rate_pix.shape[1])
            spec_rate = np.zeros(spec_shape)
            spec_variance = np.zeros(spec_shape)
            for i in np.arange(dispersion.shape[0]):
                # Background not yet completely added. Make sure there is a trace shift to be done so that we
                # don't make an expensive call to shift() if we don't have to. Use mode='nearest' to fill in new
                # pixels with background when image is shifted.
                if trace[i] != 0.0:
                    spec_rate[:, i:i + rate_pix.shape[1]] += shift(
                        electron_rate_pix[:, :, i],
                        shift=(trace[i], 0),
                        mode='nearest'
                    ) * dispersion[i]
                    spec_variance[:, i:i + rate_pix.shape[1]] += shift(
                        electron_variance_pix[:, :, i],
                        shift=(trace[i], 0),
                        mode='nearest'
                    ) * dispersion[i]
                else:
                    spec_rate[:, i:i + rate_pix.shape[1]] += electron_rate_pix[:, :, i] * dispersion[i]
                    spec_variance[:, i:i + rate_pix.shape[1]] += electron_variance_pix[:, :, i] * dispersion[i]

                # Adding background to all other pixels, unless we are asked not to.
                if add_extended_background:
                    spec_rate[:, :i] += bg_electron_rate[i] * dispersion[i]
                    spec_rate[:, i + rate_pix.shape[1]:] += bg_electron_rate[i] * dispersion[i]
                    spec_variance[:, :i] += bg_electron_variance[i] * dispersion[i]
                    spec_variance[:, i + rate_pix.shape[1]:] += bg_electron_variance[i] * dispersion[i]
        else:
            spec_shape = (rate_pix.shape[2] + rate_pix.shape[0], rate_pix.shape[1])
            spec_rate = np.zeros(spec_shape)
            spec_variance = np.zeros(spec_shape)
            for i in np.arange(dispersion.shape[0]):
                # Background not yet completely added. Make sure there is a trace shift to be done so that we
                # don't make an expensive call to shift() if we don't have to. Use mode='nearest' to fill in new
                # pixels with background when image is shifted.
                if trace[i] != 0.0:
                    spec_rate[i:i + rate_pix.shape[1], :] += shift(
                        electron_rate_pix[:, :, i],
                        shift=(0, trace[i]),
                        mode='nearest'
                    ) * dispersion[i]
                    spec_variance[i:i + rate_pix.shape[1], :] += shift(
                        electron_variance_pix[:, :, i],
                        shift=(0, trace[i]),
                        mode='nearest'
                    ) * dispersion[i]
                else:
                    spec_rate[i:i + rate_pix.shape[1], :] += electron_rate_pix[:, :, i] * dispersion[i]
                    spec_variance[i:i + rate_pix.shape[1], :] += electron_rate_pix[:, :, i] * dispersion[i]
                # Adding background to all other pixels, unless we are asked not to.
                if add_extended_background:
                    spec_rate[:i, :] += bg_electron_rate[i] * dispersion[i]
                    spec_rate[i + rate_pix.shape[1]:, :] += bg_electron_rate[i] * dispersion[i]
                    spec_variance[:i, :] += bg_electron_variance[i] * dispersion[i]
                    spec_variance[i + rate_pix.shape[1]:, :] += bg_electron_variance[i] * dispersion[i]

        # dispersion_axis determines whether wavelength is the first or second axis
        if self.dispersion_axis == 'x' or self.projection_type == 'multiorder':
            products = wave_pix_trunc, spec_rate, spec_variance
        else:
            # if dispersion is along Y, wavelength increases bottom to top, but Y index increases top to bottom.
            # flip the Y axis to account for this.
            products = wave_pix_trunc, np.flipud(spec_rate), np.flipud(spec_variance)

        return products

    def wave_eff(self, rate):
        rate_tot = np.nansum(rate, axis=0)
        a = np.sum(rate_tot * self.wave)
        b = np.sum(rate_tot)
        if b > 0.0:
            wave_eff = a / b
        else:
            wave_eff = self.wave.mean()
        wave_eff_arr = np.array([wave_eff])
        return wave_eff_arr

    def get_projection_type(self):
        return self.projection_type

    def ipc_convolve(self, rate, kernel):
        fp_pix_ipc = convolve_fft(rate, kernel, normalize_kernel=False,
                                  boundary='wrap',
                                  fftn=pyfftw.interfaces.numpy_fft.fftn,
                                  ifftn=pyfftw.interfaces.numpy_fft.ifftn)
        return fp_pix_ipc

    def get_saturation_mask(self, rate=None):
        """
        Compute a numpy array indicating pixels with hard saturation (2), soft saturation (1) and no saturation (0).

        Parameters
        ----------
        rate: None or 2D np.ndarray
            Detector plane rate image used to build saturation map from

        Returns
        -------
        mask: 2D np.ndarray
            Saturation mask image
        """
        if rate is None:
            rate = self.rate_plus_bg

        saturation_mask = np.zeros(rate.shape)

        if self.calculation_config.effects['saturation']:
            fullwell = self.det_pars['fullwell']
            exp_pars = self.current_instrument.exposure_spec
            ngroup = exp_pars.ngroup
            unsat_ngroups = exp_pars.get_unsaturated_groups(rate, fullwell)

            saturation_mask[np.where(unsat_ngroups < ngroup)] = 1
            saturation_mask[np.where(unsat_ngroups < 2)] = 2

        return saturation_mask


class CombinedSignal(object):

    """
    This class takes a set of DetectorSignal instances, combines the rates appropriately, and
    provides the other information that DetectorNoise requires.  The primary use for this is in the case
    of SOSS where the detector plane contains signals from effectively three instrument configurations, one
    for each of the visible orders of the gr700xd disperser. These signals need to be combined properly before
    being used to calculate a DetectorNoise and perform a Strategy extraction.

    WARNING: This is currently set up to only support a single aperture slice. Supporting multiple orders with
    multiple slices will require further refactoring...

    Parameters
    ----------
    signal_list: list of DetectorSignal instances
    """
    def __init__(self, signal_list):
        # each signal potentially has a different grid so we want to combine everything onto one grid. instruments
        # that combine multiple signals need to provide detector_pixels so that we know how to shift the signals
        # with respect to each other.
        self.warnings = {}
        self.detector_pixel_list = []
        self.wave_pix_list = []
        # some things are common to all signals so get them from the first one
        self.parent_signal = signal_list[0]
        self.dispersion_axis = self.parent_signal.dispersion_axis

        if len(self.parent_signal.rate_list) > 1:
            msg = "Combining multiple instrument signals onto single detector only supports a single slice."
            raise NotImplementedError(value=msg)

        maxx = 0
        maxy = 0
        for i, s in enumerate(signal_list):
            self.warnings.update(s.warnings)
            pixels = s.detector_pixels
            if pixels is None:
                inst = s.current_instrument.instrument['instrument']
                msg = "No detector pixel configuration set for instrument %s." % inst
                raise DataError(value=msg)
            self.detector_pixel_list.append(pixels)
            ny, nx = s.rate.shape
            if self.dispersion_axis == 'x':
                maxx = max(maxx, nx + pixels[0])
                maxy = max(maxy, ny)
            else:
                maxx = max(maxx, nx)
                maxy = max(maxy, ny + pixels[0])
            self.wave_pix_list.append(s.wave_pix)

        self.grid = coords.IrregularGrid(np.arange(maxy), np.arange(maxx))
        self.pixgrid_list = [self.grid]
        self.dist = self.grid.dist()
        self.det_mask = np.ones_like(self.dist)

        # these parameters are the same for each signal
        self.wave = self.parent_signal.wave
        self.total_flux = self.parent_signal.total_flux
        self.fp_rate = self.parent_signal.fp_rate
        self.bg_fp_rate = self.parent_signal.bg_fp_rate
        self.background = self.parent_signal.background
        self.flux_cube_list = self.parent_signal.flux_cube_list
        self.flux_plus_bg_list = self.parent_signal.flux_plus_bg_list

        self.aperture_list = self.parent_signal.aperture_list
        self.projection_type = self.parent_signal.projection_type
        self.current_instrument = self.parent_signal.current_instrument
        self.spatial_grid = self.parent_signal.grid
        self.det_pars = self.parent_signal.det_pars
        self.calculation_config = self.parent_signal.calculation_config
        if self.parent_signal.det_pars['rn_correlation']:
            self.read_noise_correlation_matrix = self.parent_signal.read_noise_correlation_matrix

        self.rate_list = [{
            'fp_pix': np.zeros_like(self.dist),
            'fp_pix_no_ipc': np.zeros_like(self.dist),
            'fp_pix_variance': np.zeros_like(self.dist)
        }]
        self.rate_plus_bg_list = [{
            'fp_pix': np.zeros_like(self.dist),
            'fp_pix_no_ipc': np.zeros_like(self.dist),
            'fp_pix_variance': np.zeros_like(self.dist)
        }]
        saturation = np.zeros_like(self.dist)

        for i, s in enumerate(signal_list):
            ny, nx = s.rate_list[0]['fp_pix'].shape
            # set up position within the combined rate image to add each rate
            if self.dispersion_axis == 'x':
                ly = 0
                uy = 0
                # shifts between the rates are referenced to the first rate calculated
                lx = int(self.detector_pixel_list[i][0] - self.detector_pixel_list[0][0])
                ux = int(self.rate_list[0]['fp_pix'].shape[1] - s.rate_list[0]['fp_pix'].shape[1] - lx)
            else:
                # shifts between the rates are referenced to the first rate calculated
                ly = int(self.detector_pixel_list[i][0] - self.detector_pixel_list[0][0])
                uy = int(self.rate_list[0]['fp_pix'].shape[0] - s.rate_list[0]['fp_pix'].shape[0] - ly)
                lx = 0
                ux = 0

            # this will cause pixels that don't overlap to be filled with the background. this is correct for imaging
            # and slitless, but may need revisiting for other cases.
            new_r = np.pad(s.rate_list[0]['fp_pix'], ([ly, uy], [lx, ux]), mode='edge')
            new_r_bg = np.pad(s.rate_plus_bg_list[0]['fp_pix'], ([ly, uy], [lx, ux]), mode='edge')
            new_r_noipc = np.pad(s.rate_list[0]['fp_pix_no_ipc'], ([ly, uy], [lx, ux]), mode='edge')
            new_r_bg_noipc = np.pad(s.rate_plus_bg_list[0]['fp_pix_no_ipc'], ([ly, uy], [lx, ux]), mode='edge')
            new_r_var = np.pad(s.rate_list[0]['fp_pix_variance'], ([ly, uy], [lx, ux]), mode='edge')
            new_r_bg_var = np.pad(s.rate_plus_bg_list[0]['fp_pix_variance'], ([ly, uy], [lx, ux]), mode='edge')
            new_sat = np.pad(s.saturation_list[0], ([ly, uy], [lx, ux]), mode='constant')

            saturation = np.maximum(saturation, new_sat)
            self.rate_list[0]['fp_pix'] += new_r
            self.rate_list[0]['fp_pix_no_ipc'] += new_r_noipc
            self.rate_list[0]['fp_pix_variance'] += new_r_var
            self.rate_plus_bg_list[0]['fp_pix'] += new_r_bg
            self.rate_plus_bg_list[0]['fp_pix_no_ipc'] += new_r_bg_noipc
            self.rate_plus_bg_list[0]['fp_pix_variance'] += new_r_bg_var

        self.saturation_list = [saturation]
        # SOSS only has one aperture so the on_detector rates are just the fp_pix rates.
        self.rate = self.rate_list[0]['fp_pix']
        self.rate_plus_bg = self.rate_plus_bg_list[0]['fp_pix']

    def get_saturation_mask(self, rate=None):
        """
        Compute a numpy array indicating pixels with hard saturation (2), soft saturation (1) and no saturation (0).
        This version just wraps whats implemented within DetectorSignal.

        Parameters
        ----------
        rate: None or 2D np.ndarray
            Detector plane rate image used to build saturation map from

        Returns
        -------
        mask: 2D np.ndarray
            Saturation mask image
        """
        if rate is None:
            rate = self.rate_plus_bg

        mask = self.parent_signal.get_saturation_mask(rate=rate)
        return mask

    def spectral_detector_transform(self):
        """
        Create engine API format dict section containing properties of wavelength coordinates
        at the detector plane.

        Returns
        -------
        t: dict (engine API compliant keys)
        """
        return self.parent_signal.spectral_detector_transform()

    def spectral_model_transform(self):
        """
        Create engine API format dict section containing properties of the wavelength coordinates
        used in the construction of a ConvolvedSceneCube.

        Returns
        -------
        t: dict (engine API compliant keys)
        """
        return self.parent_signal.spectral_model_transform()

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
        return self.parent_signal.cube_wcs_info()


class DetectorNoise(object):

    """
    This class provides the functionality to model the noise in a JWST observation

    Parameters
    ----------
    obs_signal: DetectorSignal or CombinedSignal instance
        The calculated signal at the detector plane
    observation: observation.Observation instance
        Configuration information for the Observation that generated ObsSignal
    """

    def __init__(self, obs_signal, observation):
        self.warnings = {}
        self.observation = observation

        self.grid = obs_signal.grid
        self.dist = obs_signal.dist
        self.det_mask = obs_signal.det_mask
        self.det_pars = obs_signal.det_pars
        self.calculation_config = obs_signal.calculation_config

        self.current_instrument = observation.instrument

        # Are we dealing with an image or spectroscopic mode?
        self.projection_type = self.current_instrument.projection_type

        # get saturation so we can calculate unsaturated groups
        if self.calculation_config.effects['saturation']:
            if 'fullwell' not in self.det_pars:
                msg = "Detector fullwell configuration missing for %s." % self.instrument.instrument
                raise DataError(value=msg)
            self.fullwell = self.det_pars['fullwell']
        else:
            self.fullwell = 1.0e99

        # if not explicitly set, assume we just need 2 groups to define a slope
        self.mingroups = self.det_pars.get("mingroups", 2)

        if self.calculation_config.noise['crs']:
            # if we are including CRs, we need a CR rate configured for the telescope. this rate is in events/s/cm^2/sr.
            if not hasattr(self.current_instrument.telescope, 'cr_rate'):
                msg = "No cosmic ray rate defined for %s." % self.current_instrument.telescope.tel_name.upper()
                raise DataError(value=msg)
            self.cr_rate = self.current_instrument.telescope.cr_rate

            # we need physical pixel size to calculate incident area in sq cm
            if 'pix_size' not in self.det_pars:
                msg = "The physical pixel size of the detector in microns must be specified to calculate the CR rate."
                raise DataError(value=msg)
            self.pix_size = self.det_pars["pix_size"]

            # need to scale the rate by number of pixels affected by each event.
            if not hasattr(self.current_instrument.telescope, 'cr_npixels'):
                msg = "No defined value for number of pixels affected per CR event for %s." % \
                    self.current_instrument.telescope.tel_name.upper()
                raise DataError(value=msg)
            self.cr_npixels = self.current_instrument.telescope.cr_npixels

            cr_rate = self.cr_rate * self.cr_npixels

            # cr_rate is in units of events/s/cm^2/sr.  convert pix_size in um to cm and multiply by 4pi steradians
            # to get cr_rate per pixel.
            self.pix_cr_rate = (self.pix_size * 1.0e-4)**2 * 4.0 * np.pi * cr_rate
        else:
            self.pix_cr_rate = 0.0

        self.var_pix_list, self.stdev_pix_list, self.rn_var_pix_list = self.basic_source_noise(obs_signal)
        self.var_pix, self.stdev_pix, self.var_rn_pix = self.on_detector()

    def on_detector(self):
        """
        Calculates the detector plane noise products

        Returns
        -------
        products: tuple
            detector_var - numpy.ndarray
                Full variance including all noise sources
            detector_stdev - numpy.ndarray
                Standard deviation (i.e. sqrt(detector_var))
            detector_rn_var - numpy.ndarray
                Variance stricly due to detector readnoise
        """
        aperture_sh = self.var_pix_list[0].shape
        n_apertures = len(self.var_pix_list)
        detector_shape = (aperture_sh[0] * n_apertures, aperture_sh[1])
        detector_var = np.zeros(detector_shape)
        detector_stdev = np.zeros(detector_shape)
        detector_rn_var = np.zeros(detector_shape)

        i = 0
        for var_pix, stdev_pix, rn_var_pix in zip(self.var_pix_list, self.stdev_pix_list, self.rn_var_pix_list):
            detector_var[i * aperture_sh[0]:(i + 1) * aperture_sh[0]:] = var_pix
            detector_stdev[i * aperture_sh[0]:(i + 1) * aperture_sh[0]:] = stdev_pix
            detector_rn_var[i * aperture_sh[0]:(i + 1) * aperture_sh[0]:] = rn_var_pix

            i += 1

        output = self.det_mask * detector_var, self.det_mask * detector_stdev, self.det_mask * detector_rn_var
        return output

    def basic_source_noise(self, obs_signal):
        """
        Calculate the noise using the full pixelated flux cube.

        Inputs
        ------
        obs_signal: ObsSignal class instance
            Class containing the detector flux plane or plane cube

        Returns
        -------
        var_pix_list: List of ndarrays
            The detector variance plane or plane cube
        stdev_pix_list: List of ndarrays
            The detector standard deviation plane or plane cube (literally the square root of var_pix_list).
        rn_var_pix_list: List of ndarrays
            Variance strictly due to detector readnoise
        """
        exp_pars = self.current_instrument.exposure_spec
        ff_electrons = self.det_pars['ff_electrons']

        var_pix_list = []
        stdev_pix_list = []
        rn_var_pix_list = []

        for rate_plus_bg in obs_signal.rate_plus_bg_list:
            slope_var_pix, slope_rn_var_pix = self.get_slope_variance(rate_plus_bg)
            rate_per_pix = rate_plus_bg['fp_pix']
            """
            The flat field error is a division by ~1 (the flat field is normalized), with a variance of 1/ff_electrons.
            Note that the value of the flat field response is constant for multiple ramps and multiple integrations, so
            nramps > 1 does not decrease the residual flat field noise. Due to that, a user will either have to improve the
            flat field or dither with > 1 pixel offsets. The most apparent effect for everyday ETC use is that this sets an
            upper limit on the achievable signal-to-noise ratio.

            The pixel variance upon division with a normalized flat field constructed with ff_electrons (it is assumed that
            the flat field is ideal):

            s^2(R/FF) = s^2(R) + R^2/FF_electrons
            """
            var_pix = slope_var_pix / exp_pars.nramps
            rn_var_pix = slope_rn_var_pix / exp_pars.nramps

            # Add the flat field residual noise if requested
            if self.calculation_config.noise['ffnoise']:
                var_pix += rate_per_pix ** 2 / ff_electrons

            stdev_pix = np.sqrt(var_pix)

            var_pix_list.append(var_pix)
            stdev_pix_list.append(stdev_pix)
            rn_var_pix_list.append(rn_var_pix)

        products = var_pix_list, stdev_pix_list, rn_var_pix_list
        return products

    def calc_cr_loss(self, ngroups):
        """
        Calculate the effective loss of exposure time due to cosmic rays. This uses the cosmic ray (CR) event rate
        that is contained within the telescope configuration to calculate the odds of a cosmic rays
        hitting a pixel. Pixels that are hit by cosmic rays are assumed to have ramps that are valid before the
        CR, but not after. The exposure specification is used to adjust the expectation value of the exposure
        time (i.e. ngroups) for a CR-truncated ramp.

        See discussion and links in:

        - https://confluence.stsci.edu/pages/viewpage.action?spaceKey=JWST&title=2014-11-10+Effect+of+CRs+on+SNR
        - Robberto (2010) Technical Report JWST-STScI-001928

        for more information on implementation and the numbers used.

        Parameters
        ----------
        ngroups: int
            Number of groups over which to calculate CR losses.  This will usually be number of unsaturated groups.

        Returns
        -------
        cr_ngroups: float
            This is the input ngroups scaled by the mean loss of time due to cosmic ray events.
        """
        if not self.calculation_config.noise['crs']:
            # if we're not correcting for CRs, just pass ngroups back
            cr_ngroups = ngroups
        else:
            exp_pars = self.current_instrument.exposure_spec

            # this is the average fraction of the ramp that is lost upon a CR event
            if ngroups < self.mingroups:
                # if for some reason (e.g. using single read noise model) ngroup is less than the minimum needed
                # to get a good ramp fit, then the whole ramp is lost upon a CR event
                ramp_frac = 1.0
            else:
                # the more groups you have per ramp, the less of the ramp you lose to CRs on average.
                # this formalism takes into account that there's always a fraction that's totally lost.
                # in the limit of infinite reads this converges on half the ramp since CRs are evenly
                # distributed.
                ramp_frac = 1.0 - 0.5 * (ngroups - self.mingroups) / (ngroups - 1.0)

            # the effective ramp exposure time, t_eff, is t_tot * (1 - ramp_frac * pix_cr_rate * t_ramp).
            # we recalculate t_ramp using the input ngroups to reflect that we might be using one truncated by saturation.
            t_ramp = exp_pars.tframe * ((ngroups * exp_pars.nframe) +
                     (ngroups - 1) * exp_pars.nskip + exp_pars.nextra)
            cr_ngroups = (1.0 - ramp_frac * self.pix_cr_rate * t_ramp) * ngroups

        return cr_ngroups

    def get_slope_variance(self, rate):
        """
        Compute the slope variance, based on detector properties and
        exposure specifications.

        Inputs
        ------
        rate_per_pix: ndarray
           The measured slope of the ramp (ie, the rate) per pixel
           This can be a one-dimensional wavelength spectrum (for a 1D ETC)
        or a three-dimensional cube ([wave,x,y] for a 3D ETC).

        Returns
        -------
        slope_var: ndarray
           The associated variance of the measured slope
        slope_rn_var: ndarray
           The associated variance of the readnoise only

        """
        dark_current = self.det_pars.get('dark_current', 0.0)
        rn = self.det_pars.get('rn', 0.0)
        var_fudge = self.det_pars.get('var_fudge', 1.0)
        rn_fudge = self.det_pars.get('rn_fudge', 1.0)

        exp_pars = self.current_instrument.exposure_spec

        unsat_ngroups = exp_pars.get_unsaturated_groups(rate['fp_pix_no_ipc'], self.fullwell, hard_saturation=self.mingroups)
        # scale the unsaturated ngroups by the CR loss.
        if self.calculation_config.noise['crs']:
            unsat_ngroups = np.vectorize(self.calc_cr_loss)(unsat_ngroups)

        slope_var, slope_rn_var = exp_pars.slope_variance(rate, dark_current, rn, unsat_ngroups,
                                                          var_fudge, rn_fudge)

        return slope_var, slope_rn_var

    def get_readnoise_slope_variance(self, rate):
        """
        Compute the variance of the read noise only, excluding all other noise contributions.

        Inputs
        ------
        none

        Returns
        -------
        var_rn: ndarray
            The variance of the readnoise given the current detector slope parameters. Note that
            if the readnoise is switched off, the associated variance will of course be zero.
        """
        if self.calculation_config.noise['readnoise']:
            rn = self.det_pars['rn']
        else:
            rn = 0.0

        exp_pars = self.current_instrument.exposure_spec
        unsat_ngroups = exp_pars.get_unsaturated_groups(rate['fp_pix_no_ipc'], self.fullwell, hard_saturation=self.mingroups)
        # scale the unsaturated ngroups by the CR loss.
        if self.calculation_config.noise['crs']:
            unsat_ngroups = np.vectorize(self.calc_cr_loss)(unsat_ngroups)

        var_rn = exp_pars.rn_variance(rn, unsat_ngroups=unsat_ngroups)

        return var_rn


def calculate_sn(inp, webapp=False):
    """
    This is a function to do the 'forward' exposure time calculation where given a dict
    in engine API input format we calculate the resulting Signal/Noise and return a Report
    on the results.

    Parameters
    ----------
    input: dict
        Engine API format dictionary containing the information required to perform the calculation.
    psf_ibrary : psf_library.PSFLibrary instance
        Library of PSF files (e.g. produced by webbpsf) to be used in the calculation

    Returns
    -------
    report.Report instance
    """
    warnings = {}
    try:
        scene_configuration = inp['scene']
        background = inp['background']
        instrument_configuration = inp['configuration']
        strategy_configuration = inp['strategy']
    except KeyError as e:
        message = "Missing information required for the calculation: %s" % str(e)
        raise EngineInputError(value=message)

    # get the calculation configuration from the input or use the defaults
    if 'calculation' in inp:
        calc_config = CalculationConfig(config=inp['calculation'])
    else:
        calc_config = CalculationConfig()

    # #### BEGIN calculation #### #
    """
    This section currently implements the Pandeia engine API.  As the engine's object model
    is refactored, this section will have to change accordingly.
    """
    # check for empty scene configuration and set it up properly if it is empty.
    if len(scene_configuration) == 0:
        scene_configuration = build_empty_scene()

    scene = Scene(input=scene_configuration, webapp=webapp)
    warnings.update(scene.warnings)
    instrument = InstrumentFactory(config=instrument_configuration, webapp=webapp)
    warnings.update(instrument.warnings)
    strategy = StrategyFactory(instrument, config=strategy_configuration, webapp=webapp)

    # strategies can have different figures of interest that need to be calculated.
    # in most cases, S/N is what is desired.  however, for coronagraphy the figure of interest
    # is sometimes the contrast that can be achieved.  in this case, strategy.calc_type will be
    # set to 'contrast' and we need to run calculate_contrast.
    if hasattr(strategy, "calc_type"):
        if strategy.calc_type == "contrast":
            r = calculate_contrast(inp, webapp=webapp)
        else:
            msg = "Unsupported calculation type: %s" % strategy.calc_type
            raise EngineInputError(value=msg)
    else:
        # set up the observation and then do S/N calculation...
        obs = observation.Observation(
            scene=scene,
            instrument=instrument,
            strategy=strategy,
            background=background,
            webapp=webapp
        )

        # seed the random number generator
        seed = obs.get_random_seed()
        np.random.seed(seed=seed)

        # Sometimes there is more than one exposure involved so implement lists for signal and noise
        my_detector_signal_list = []
        my_detector_noise_list = []
        my_detector_saturation_list = []

        if hasattr(strategy, 'dithers'):
            dither_list = strategy.dithers
        else:
            dither_list = [{'x': 0.0, 'y': 0.0}]

        for dither in dither_list:
            o = copy.deepcopy(obs)
            o.scene.offset(dither)
            # Calculate the signal rate in the detector plane. If they're configured, need to loop through
            # configured orders to include all dispersed signal.
            if instrument.projection_type == 'multiorder':
                norders = instrument.disperser_config[instrument.instrument['disperser']]['norders']
                orders = list(range(1, norders + 1))
            else:
                orders = None

            if orders is not None:
                order_signals = []
                for order in orders:
                    order_signals.append(DetectorSignal(o, calc_config=calc_config, webapp=webapp, order=order))
                my_detector_signal = CombinedSignal(order_signals)
            else:
                my_detector_signal = DetectorSignal(o, calc_config=calc_config, webapp=webapp, order=None)

            my_detector_noise = DetectorNoise(my_detector_signal, o)

            # Every dither has a saturation map
            my_detector_saturation = my_detector_signal.get_saturation_mask()

            my_detector_signal_list.append(my_detector_signal)
            my_detector_noise_list.append(my_detector_noise)
            my_detector_saturation_list.append(my_detector_saturation)

        # Use the strategy to get the extracted signal/noise products
        extracted_sn = strategy.extract(my_detector_signal_list, my_detector_noise_list)
        warnings.update(extracted_sn['warnings'])
        # #### END calculation #### #
        r = Report(inp, my_detector_signal_list, my_detector_noise_list, my_detector_saturation_list, extracted_sn, warnings)

    return r


def calculate_contrast(inp, webapp=False):
    """
    This is a function to do the 'forward' exposure time calculation where given a dict
    in engine API input format we calculate the resulting coronagraphic contrast and return a Report
    on the results.

    While this method is meant for coronagraphic modes, it will work also for regular imaging modes.

    Parameters
    ----------
    input: dict
        Engine API format dictionary containing the information required to perform the calculation.
    psf_ibrary : psf_library.PSFLibrary instance
        Library of PSF files (e.g. produced by webbpsf) to be used in the calculation

    Returns
    -------
    report.Report instance
    """
    warnings = {}
    try:
        scene_configuration = inp['scene']
        background = inp['background']
        instrument_configuration = inp['configuration']
        strategy_configuration = inp['strategy']
        psf_subtraction_configuration = inp['strategy']['psf_subtraction_source']
    except KeyError as e:
        message = "Missing information required for the calculation: %s" % str(e)
        raise EngineInputError(value=message)

    # get the calculation configuration from the input or use the defaults
    if 'calculation' in input:
        calc_config = CalculationConfig(config=inp['calculation'])
    else:
        calc_config = CalculationConfig()

    # #### BEGIN calculation #### #
    """
    This section implements the Pandeia engine API.
    """
    # check for empty scene configuration and set it up properly if it is empty.
    if len(scene_configuration) == 0:
        scene_configuration = build_empty_scene()

    instrument = InstrumentFactory(config=instrument_configuration, webapp=webapp)
    warnings.update(instrument.warnings)
    # move the psf_reference to a pre-determined and fixed location outside of the FOV
    strategy = StrategyFactory(instrument, config=strategy_configuration, webapp=webapp)

    # Check for user-specified dithers (the contrast calculation will add an additional 2 fictional dithers).
    if not hasattr(strategy, 'dithers'):
        if len(strategy['dithers'] != 1):
            message = "Contrast calculations currently only support a single dither to be passed in the strategy: %s" % str(e)
            raise EngineInputError(value=message)

    psf_subtraction_xy = strategy.psf_subtraction_xy
    pointing_error = strategy.pointing_error
    psf_subtraction_configuration['position']['x_offset'] = psf_subtraction_xy[0]
    psf_subtraction_configuration['position']['y_offset'] = psf_subtraction_xy[1]

    # add the psf_reference to the scene
    scene_configuration.append(psf_subtraction_configuration)
    scene = Scene(input=scene_configuration, webapp=webapp)
    if hasattr(strategy, "scene_rotation"):
        scene.rotate(strategy.scene_rotation)
    warnings.update(scene.warnings)

    # Add the appropriate dithers for the PSF reference star and the unocculted dither
    psf_subtraction_dither = {
        'x': -psf_subtraction_xy[0] - pointing_error[0],
        'y': -psf_subtraction_xy[1] - pointing_error[1]
    }
    unocculted_dither = {'x': strategy.unocculted_xy[0], 'y': strategy.unocculted_xy[1]}

    strategy.dithers.append(psf_subtraction_dither)
    strategy.dithers.append(unocculted_dither)
    strategy.on_target = [True, False, False]

    # set up the observation...
    obs = observation.Observation(
        scene=scene,
        instrument=instrument,
        strategy=strategy,
        background=background,
        webapp=webapp
    )

    # seed the random number generator
    seed = obs.get_random_seed()
    np.random.seed(seed=seed)

    # Sometimes there is more than one exposure involved so implement lists for signal and noise
    my_detector_signal_list = []
    my_detector_noise_list = []
    my_detector_saturation_list = []

    for dither in strategy.dithers:
        # make a new deep copy of the observation for each dither so that each position is offset
        # from the center position. otherwise the offsets get applied cumulatively via the reference.
        o = copy.deepcopy(obs)
        o.scene.offset(dither)
        # Calculate the signal rate in the detector plane
        my_detector_signal = DetectorSignal(o, calc_config=calc_config, webapp=webapp)
        my_detector_noise = DetectorNoise(my_detector_signal, o)

        # Every dither has a saturation map
        my_detector_saturation = my_detector_signal.get_saturation_mask()

        my_detector_signal_list.append(my_detector_signal)
        my_detector_noise_list.append(my_detector_noise)
        my_detector_saturation_list.append(my_detector_saturation)

    # We need a regular S/N of the target source
    extracted_sn = strategy.extract(my_detector_signal_list, my_detector_noise_list)
    warnings.update(extracted_sn['warnings'])

    # Use the strategy to get the extracted contrast products
    grid = my_detector_signal_list[0].grid

    aperture = strategy.aperture_size
    annulus = strategy.sky_annulus

    # Create a list of contrast separations for which to calculate the contrast
    bounds = grid.bounds()
    ncontrast = strategy.ncontrast
    contrasts = np.zeros(ncontrast)
    contrast_separations = np.linspace(0 + aperture, bounds['xmax'] - annulus[1], ncontrast)
    contrast_azimuth = np.radians(strategy.contrast_azimuth)
    contrast_xys = [(separation * np.sin(contrast_azimuth),
                     separation * np.cos(contrast_azimuth)) for separation in contrast_separations]

    # Calculate contrast at each separation
    for i, contrast_xy in enumerate(contrast_xys):
        strategy.target_xy = contrast_xy
        extracted = strategy.extract(my_detector_signal_list, my_detector_noise_list)
        contrasts[i] = extracted['extracted_noise']

    # What is the flux of the unocculted star.
    # We set the dither weights such that only the unocculted dither is used.
    strategy.dither_weights = [0, 0, 1]
    strategy.on_target = [False, False, True]
    strategy.target_xy = strategy.unocculted_xy
    extract_unocculted = strategy.extract(my_detector_signal_list, my_detector_noise_list)

    # when a source is offset to unocculted_xy, it can be bright enough to cause saturation
    # flags to be raised.  however, since this is an "artifactual" offset, those saturation
    # flags are bogus. the hackish fix is to pop this bogus saturation map off the list
    # and append a new one filled with zeros.
    bogus_sat = my_detector_saturation_list.pop()
    my_detector_saturation_list.append(np.zeros(bogus_sat.shape))

    # Contrast is relative to the unocculted on-axis star.
    contrasts /= extract_unocculted['extracted_flux']
    contrast_curve = [contrast_separations, contrasts]

    # #### END calculation #### #

    # Add the contrast curve and link relevant saturation maps to the extracted_sn dict for passing
    extracted_sn['contrast_curve'] = contrast_curve

    r = Report(inp, my_detector_signal_list, my_detector_noise_list, my_detector_saturation_list, extracted_sn, warnings)
    return r


def calculate_exposure_time(inp, webapp=False, **kwargs):
    """
    This is a function to do the 'reverse' exposure time calculation where given a desired
    signal-to-noise ratio we calculate the optimal exposureSpecification and return a Report
    on the results.

    Parameters
    ----------
    input: dict
        Dictionary containing the information required to perform the calculation.
    PDFLibrary : psf_library.PSFLibrary instance
        Library of PSF files (e.g. produced by webbpsf) to be used in the calculation

    Returns
    -------
    report.Report instance
    """
    raise NotImplementedError("Reverse ETC calculations are not yet implemented.")
