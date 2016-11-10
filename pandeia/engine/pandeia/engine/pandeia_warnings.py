# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import


"""
Warning defined here are for conditions that can affect the results of a calculation,
but do not prevent a calculation from proceeding without error.

API warnings about the use of default values in calculations are handled programmatically
within each class.
"""

# these are warning messages that are common to all modules. the specific warnings currently
# defined will probably go away once API and sanity checking is fully implemented and verified
# throughout the system. these messages get combined with each of the module-specific messages using
# the update() method.
standard_warning_messages = {
    "no_sanity_check": "No sanity checking implemented for %s: %s",
    "no_api_check": "No API checking implemented for %s: %s"
}

# warning messages specific to Strategy and its sub-classes.
strategy_warning_messages = {
    "background_region_too_small": "Background region smaller than source extraction region. This can adversely affect the SNR.",
    "extraction_aperture_truncated": "Extraction aperture partially outside of the field of view.",
    "background_region_truncated": "Background estimation region partially outside of the field of view.",
    "extraction_aperture_undersampled": "Extraction aperture undersampled. Pixel area %.2f%% less than requested area.",
    "background_region_undersampled": "Background estimation region undersampled. Pixel area %.2f%% less than requested area.",
    "upper_background_region_missing": "Upper background estimation region outside of field of view.",
    "upper_background_region_truncated": "Upper background estimation region truncated by %.2f%%.",
    "lower_background_region_missing": "Lower background estimation region outside of field of view.",
    "lower_background_region_truncated": "Lower background estimation region truncated by %.2f%%.",
    "target_x_position_not_used": "Target X position is not meaningful and thus not used in spectral aperture photometry.",
    "target_occulted": "Specified target position is occulted by the coronagraphy mask."
}
strategy_warning_messages.update(standard_warning_messages)

# warning messages specific to Normalization and its sub-classes
normalization_warning_messages = {
    "normalized_to_zero_flux": "Zero flux at reference wavelength. Spectrum left unscaled.",
    "unsupported_normalization_bandpass": "Bandpass specification, %s, not currently supported, but may work."
}
normalization_warning_messages.update(standard_warning_messages)

# warning messages specific to SED and its sub-classes. currently none defined...
sed_warning_messages = {}
sed_warning_messages.update(standard_warning_messages)

# warning messages specific to Telescope and its sub-classes
telescope_warning_messages = {
    "telescope_ote_efficiency_missing": "Telescope OTE throughput mis-configured or unavailable. Using default value of 1.0.",
    "telescope_background_missing": "Telescope notional background mis-configured or unavailable. Using default value of 0.0."
}
telescope_warning_messages.update(standard_warning_messages)

# warning messages specific to AstroSpectrum and its classes and methods
astrospectrum_warning_messages = {
    "max_scene_size_reached": "Scene requires a total field-of-view of %.3f arcsec. Using the configured maximum of %.3f arcsec.",
    "scene_fov_too_small": "Field-of-view size, %.3f, too small to encompass any of the defined sources.",
    "wavelength_truncated_blue": "Spectrum blue limit, %.2f, does not extend to instrument configuration blue limit, %.2f.",
    "wavelength_truncated_red": "Spectrum red limit, %.2f, does not extend to instrument configuration red limit, %.2f.",
    "scene_range_truncated": "Combined wavelength range of scene [%.2f, %.2f] less than instrument configuration's range [%.2f, %.2f].",
    "spectrum_missing_blue": "Spectrum [%.2f, %.2f] does not extend to instrument configuration blue limit, %.2f.",
    "spectrum_missing_red": "Spectrum [%.2f, %.2f] does not extend to instrument configuration red limit, %.2f."
}
astrospectrum_warning_messages.update(standard_warning_messages)
