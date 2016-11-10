# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import numpy as np
import numpy.ma as ma


class ExposureSpecification:

    """
    Encapsulates all data associated with exposure specification
    """

    def __init__(self, pattern, ngroup, nint, nexp, tframe,
                 nframe=1, subarray='Full', nskip=0, nextra=0):
        """
        Create a generic Exposure Specification.

        The equation for on-source time comes from Appendix D of the JWST PPS Design
        proposal instructions at:

        http://www.stsci.edu/institute/org/oed/JWST_PPS_Design/proposal_instructions_folder/AppD-Photontime.doc

        Inputs
        ------
        pattern: str
            name of readout pattern

        ngroup: int
            number of groups per integration

        nint: int
            number of integrations per exposure

        nexp: int
            number of exposures

        tframe: float
            number of seconds per frame

        nframe: int (optional, default 1)
            number of frames per group

        subarray: str (optional, default 'Full')
            name of subarray to be read.
            If None, readout is full frame

        nskip: int (optional, default 0)
            number of skipped frames per group. Only
            supported by some readout patterns.

        nextra: int (optional, default 0)
            number of extra frames read out per integration. as of now, this will
            be 1 in the case of NIRISS and FGS full frame exposures and 0 in all
            other cases.
        """

        # Eventually We might do some value checking at this point
        # to ensure that the specified values are consistent with the
        # named pattern.

        # Required parameters
        self.pattern = pattern
        self.ngroup = ngroup
        self.nint = nint
        self.nexp = nexp
        self.tframe = tframe

        # Optional parameters
        self.nframe = nframe
        self.subarray = subarray
        self.nskip = nskip
        self.nextra = nextra

        # Derived quantities
        self.nramps = self.nint * self.nexp
        self.tgroup = self.tframe * (self.nframe + self.nskip)
        self.tramp = self.tframe * ((self.ngroup * self.nframe) +
                                    (self.ngroup - 1) * self.nskip + self.nextra)
        self.on_source_time = self.tramp * self.nexp * self.nint

    def get_unsaturated_groups(self, slope, fullwell, hard_saturation=2):
        """
        Calculate the number of unsaturated groups in each pixel, given a specified slope and full well value.
        The formula is calculated by equating the frame time to the full well value and isolating ngroup. The
        resulting fractional value is then rounded down to an integer. There is a minimum sensible value of
        unsaturated number of groups that defines hard saturation. The default is 2, but this can potentially
        be set to a higher value, or even 1(!) for certain observing modes. It is not expected that a value of
        0 can ever be anything but saturated.

        Parameters
        ----------
        slope: ndarray
            The measured slope of the ramp (ie, the rate) per pixel
            This can be a one-dimensional wavelength spectrum (for a 1D ETC)
            or a three-dimensional cube ([wave,x,y] for a 3D ETC).

        fullwell: positive integer
            The number of electrons defining a full well (beyond which the pixel reads become unusuable for science
            due to saturation).

        hard_saturation: positive integer
            The minimum number of groups allowed to define an unsaturated measurement.

        Returns
        -------
        unsat_ngroups: MaskedArray
            The number of unsaturated groups present in the ramp. The mask separates pixels that have hard saturation
            from those that do not. The former will not have a valid noise measurement and cannot be used in a strategy.

        """

        """
        We add a very small number to the slope here avoid dividing by zero. Before pixels with 0 slope were masked,
        but that is incorrect. 0 slope is still a valid number.
        """
        unsat_ngroups_frac = ((fullwell / (self.tframe * (slope + 1e-10))) + self.nskip -
                              self.nextra) / (self.nframe + self.nskip)

        max_ngroups = np.floor(unsat_ngroups_frac)
        max_ngroups = max_ngroups.clip(0, self.ngroup)

        unsat_ngroups = ma.masked_less(max_ngroups, hard_saturation)

        return unsat_ngroups

    def slope_variance(self, rate, dark_current, readnoise, unsat_ngroups,
                       var_fudge=1.0, rn_fudge=1.0):
        """
        Calculate the variance of a specified MULTIACCUM slope.

        Inputs
        ------
        rate: ndarray
            The measured slope of the ramp (ie, the rate) per pixel
            This can be a one-dimensional wavelength spectrum (for a 1D ETC)
            or a three-dimensional cube ([wave,x,y] for a 3D ETC).

        dark_current: float
            Dark current (electrons/s).

        readnoise: float
            Readnoise per pixel.

        unsat_ngroups: ndarray
            Number of unsaturated groups for each pixel

        var_fudge: float
            Fudge factor to apply to variance to match IDT results

        rn_fudge: float
            Fudge factor to apply to readnoise to match IDT results

        Returns
        -------
        slope_var: ndarray
            Variance associated with the input slope.
        slope_rn_var: ndarray
            The associated variance of the readnoise only
        """

        # Rename variables for ease of comparison with Robberto (35).
        # The noise calculation depends on the pixel rate BEFORE IPC convolution. fp_pix_variance takes that pre-IPC rate
        # and scales it by the quantum yield and Fano factor to get the per-pixel variance in the electron rate.
        variance_per_pix = rate['fp_pix_variance']

        rn = readnoise
        n = unsat_ngroups  # we discard any saturated groups
        m = self.nframe + self.nskip
        tframe = self.tframe
        tgroup = self.tgroup

        # Compute the variance of a MULTIACCUM slope using Robberto's formula (35). The rn_variance is also
        # calculated using the unsaturated number of groups.
        slope_rn_var = self.rn_variance(rn, unsat_ngroups=n, rn_fudge=rn_fudge)

        # The final slope variance slope may also be worse than the theoretical best value.
        # (also from Glasse et al. 2015, PASP 127 686).
        if var_fudge != 1:
            slope_var = (6. / 5.) * (n ** 2. + 1.) / (n * (n ** 2. - 1.)) * \
                ((variance_per_pix * var_fudge + dark_current) / tgroup) * \
                (1. - (5. / 3.) * (m ** 2. - 1.) / (m * (n ** 2. + 1.)) * (tframe / tgroup)) + \
                slope_rn_var
        else:
            slope_var = (6. / 5.) * (n ** 2. + 1.) / (n * (n ** 2. - 1.)) * \
                ((variance_per_pix + dark_current) / tgroup) * \
                (1. - (5. / 3.) * (m ** 2. - 1.) / (m * (n ** 2. + 1.)) * (tframe / tgroup)) + \
                slope_rn_var

        # The default fill value for masked arrays is a finite number, so convert to ndarrays, and fill with NaNs to make
        # sure missing values are interpreted as truly undefined downstream.

        slope_var = ma.filled(slope_var, fill_value=np.nan)
        slope_rn_var = ma.filled(slope_rn_var, fill_value=np.nan)

        return slope_var, slope_rn_var

    def rn_variance(self, readnoise, unsat_ngroups=None, rn_fudge=1.0):
        """
        Calculate the variance due to read noise only.

        Inputs
        ------
        unsat_ngroups: ndarray
           The number of unsaturated groups.

           If unsat_ngroups not supplied, the approximation is that the read noise is
           negligible for pixels with signal rates high enough to saturate in part of the ramp.
           This is probably always a good approximation.

        readnoise: float
          Readnoise per pixel.

        Returns
        -------
        var_rn: ndarray if unsat_ngroups is supplied, float if unsat_ngroups is not supplied.
           Variance associated with the read noise only.

        """
        if unsat_ngroups is None:
            n = self.ngroup
        else:
            n = unsat_ngroups

        rn = readnoise
        m = self.nframe + self.nskip
        tgroup = self.tgroup

        var_rn = 12. * rn ** 2. / (m * n * (n ** 2. - 1.) * tgroup ** 2.)

        # The readnoise on the slope may be worse than the theoretical best value
        # (see Glasse et al. 2015, PASP 127 686).
        if rn_fudge != 1:
            var_rn *= rn_fudge

        return var_rn
