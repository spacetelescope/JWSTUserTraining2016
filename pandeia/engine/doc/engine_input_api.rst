Currently Supported Calculation Modes (03-2016)
===============================================

The following following JWST instrument/mode combinations are currently implemented, working,
and part of the nightly regression tests:

* ``miri:``
    - imaging
    - mrs
    - coronagraphy
    - lrsslit
    - lrsslitless

* ``nircam:``
    - sw_imaging
    - lw_imaging
    - ssgrism
    - wfgrism
    - coronagraphy

* ``nirspec:``
    - ifu
    - msa
    - fixed_slit

* ``niriss:``
    - imaging
    - soss
    - ami
    - wfss

The following HST instrument/mode combination is not yet fully implemented, but will be available
for engine testing and validation:

* ``wfc3:``
    - imaging

The following WFIRST instrument/mode combination is now implemented with very notional reference data:

* ``wfirstimager:``
    - imaging
    - ifu

Overview of Inputs
==================

The engine input api is a dict.  At the top level of the dict, there are:

background:
    string for notional background lookup, or see below

configuration:
    dict: camera/telescope configuration

scene:
    list: list of sources

strategy:
    dict: strategy configuration


calculation:
    dict: test interface

fake_exception:
    test interface

error, server_test:
    reserved by server


Details of Inputs for ETC Calculation
=====================================

NOTE - keep all of the following in synch with: pandeia/engine/helpers/schema

As of now, the engine requires the following information to perform a
calculation:

scene: list (no default)
  This is a list of Sources. Each Source is described by a dict with
  the following keys:

    position: dict
      Source position parameters described by the following keys:

        x_offset: float (default 0.0)
            Detector plane X offset from FOV center in arcseconds. Positive to the right.
        y_offset: float (default 0.0)
            Detector plane Y offset from FOV center in arcseconds. Positive is up.
        orientation: float (default 0.0)
            Detector plane orientation in degrees. Positive is in direction of +X.
            (e.g. orientation=90 is UP, and orientation=-90 is DOWN)

    shape: dict
      Source shape parameters described by the following keys:

        geometry: string (default "point")
            Supported geometries are "point", "gaussian2d", "flat", and "sersic". Required
            additional source shape parameters are contingent on this parameter:

                "point" requires no additional parameters.

                "gaussian2d", "flat", and "sersic" all require these parameters:
                    major: float (default 0.1)
                        Semi-major axis in arcseconds. For "flat" this sets the size, for "gaussian2d"
                        this sets the sigma, and for "sersic" this sets a scale length where
                        I(r) = I(0)/e.
                    minor: float (default 0.1)
                        Semi-minor axis in arcseconds

                "sersic" requires one additional parameter:
                    sersic_index: float (default 1.0)
                        Power law index that sets the shape of a sersic profile.
                        sersic_index = 1.0 --> exponential
                        sersic_index = 0.5 --> gaussian
                        sersic_index = 4.0 --> de Vaucouleurs

    spectrum: dict
      Source spectral parameters described by the following keys:

        redshift: float (default 0.0)
            Redshift to apply to the continuum. Since lines are added with physical units for their strength,
            they are added to the spectrum after normalization and redshift.

        extinction: dict
          Defines how the spectrum is reddened by interstellar dust

            law: string
                Extinction law to use. Supported laws are
                    * ``mw_rv_31`` - WD01 Milky Way curve for an R_V value of 3.1 (default)
                    * ``mw_rv_40`` - WD01 Milky Way curve for an R_V value of 4.0
                    * ``mw_rv_55`` - WD01 Milky Way curve for an R_V value of 5.5
                    * ``hd210121`` - WD01 Extinction curve for high-latitude molecular cloud hd210121 with C/H = b_C = 40 ppm
                                     in log-normal size dists
                    * ``lmc_avg``  - WD01 Average extinction curve for the LMC with C/H = b_C = 20 ppm in log-normal size dists
                    * ``lmc_2``    - WD01 LMC extinction curve with C/H = b_C = 10 ppm in log-normal size dists (30 Dor region)
                    * ``smc_bar``  - WD01 Extinction curve in SMC bar with C/H = b_C = 0 ppm in log-normal size dists
                    * ``chapman09`` - Chapman et al. (2009) mid-IR extinction curve derived from three molecular clouds:
                                      Ophiuchus, Perseus, and Serpens
            value: float
                Level of extinction in units of unit
            unit: string
                Units of extinction.  Allowed values are ``nh`` for hydrogen column density (cm^-2) and "mag" for magnitudes
                of extinction in specified bandpass, ext_bandpass
            bandpass: string
                Bandpass to which extinction is normalized to if unit="mag".  Allowed values are v, j, h, and k.

        normalization: dict
          Defines how the spectrum is to be scaled.

            type: string
                Method of normalization to perform.  Supported methods are
                    * ``at_lambda`` - Specify norm_flux in fluxunit at a specfic wavelength, norm_wave
                    * ``hst`` - Specify a bandpass in the form of an "obsmode" string to pass along to pysynphot along with fluxunit and norm_flux
                    * ``jwst`` - Specify a bandpass as an instrument configuration in the form of a comma-separated string <instrument>,<mode>,<filter> along with fluxunit and norm_flux
                    * ``photsys`` - Specify bandpass in the form of a comma-separated string <photsys>,<filter>
                    * ``none`` - Do not normalize spectrum.  Only valid for a spectrum type of 'input'.

            norm_wave: float
                Reference wavelength in 'norm_waveunit' at which spectrum will be scaled for type 'at_lambda'.
                Ignored for other normalization types.
            norm_flux: float
                Reference flux in 'norm_fluxunit' to which spectrum will be scaled.
            norm_fluxunit: string
                Specify the flux units in which the normalization should occur.
                Supports flam, fnu, vegamag, abmag, mjy, ujy, njy, jy
            norm_waveunit: string
                Specify the wavelength units used in normalization
            bandpass: string
                Specifies the key used to obtain the normalization bandpass for
                types 'hst', 'jwst', and 'photsys'.

        sed: dict
          Defines the spectral energy distribution of the spectrum.

            sed_type: string
                Type of the spectral energy distribution. Each type requires its own set
                of parameters. The analytic sed_type's (none, flat, powerlaw, flat) all
                require 'wmin', 'wmax', and 'sampling' to define the range and wavelength
                sampling over which the model spectrum is calculated. However, they are only
                available in the API for testing purposes and should not be configured via
                the UI.

                    **no_continuum** - No continuum, specifically Flux = 0.0 over specified range [wmin, wmax]
                        wmin: float (default 0.5)
                            Minimum wavelength in microns
                        wmax: float (default 30.0)
                            Maximum wavelength in microns
                        sampling: int (default 200)
                            Sets the logarithmic wavelength sampling of the model spectrum

                    **flat** - Flat spectrum in specified units calculated over specified range [wmin, wmax]
                        wmin: float (default 0.5)
                            Minimum wavelength in microns
                        wmax: float (default 30.0)
                            Maximum wavelength in microns
                        sampling: int (default 200)
                            Sets the logarithmic wavelength sampling of the model spectrum
                        unit: string
                            Units of spectrum, either 'fnu' or 'flam'

                    **powerlaw** - Powerlaw spectrum where F ~ lambda ^ index calculated over range [wmin, wmax]
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

                    **blackbody** - Blackbody spectrym calculated over range [wmin, wmax]
                        wmin: float (default 0.5)
                            Minimum wavelength in microns
                        wmax: float (default 30.0)
                            Maximum wavelength in microns
                        sampling: int (default 200)
                            Sets the logarithmic wavelength sampling of the model spectrum
                        temp: float
                            Temperature of the blackbody

                    **phoenix** - Parameterized stellar atmosphere models calculated by the Phoenix group
                        key: string
                            In webapp mode, a key is used to look up a predefined set of parameters. If not
                            in webapp mode and if key is not provided, model parameters can be passed directly:
                        teff: float
                            Effective temperature. Allowed range is 2000 K to 70000 K
                        log_g: float
                            Surface gravity in log10(cgs) units. Allowed range is 0.0 to 5.5.
                        metallicity: float
                            Metallicity in units of log10(solar metallicity). Allowed range is -4.0 to 0.5.

                    **hst_calspec** - HST standard star spectra
                        key: string
                            Key used to look up which spectrum load.

                    **galaxies** - Integrated spectra of galaxies from Brown et al. (2014)
                        key: string
                            Key used to look up which spectrum load.

                    **input** - spectrum provided via input arrays
                        spectrum: list-like or numpy.ndarray
                            The 0th index is taken to be wavelength in units of 'mJy'.
                            The 1st index is taken to be the flux in units of 'microns'.

        lines: list (default [])
          List of line definitions. Each definition is a dict with keys:

              name: string (default 'no name')
                  Name of line (e.g. 'Hydrogen Alpha')
              center: float (default 5.2)
                  Wavelength at line center in w_unit
              strength: float (default 1.0e-14)
                  Strength of line in erg/cm^2/s for emission or
                  optical depth for absorption
              profile: string
                  Line profile type:
                    * gaussian      *default*
                    * voigt          NOT YET IMPLEMENTED
              emission_or_absorption: string
                  Line type:
                    * emission      *default*
                    * absorption

            A profile type of **gaussian** requires one additional parameter:

              width: float (default 200.0)
                  Full-width half-max of line in km/s

            When implemented, profile type of **voigt** will require two additional parameters:

              gaussian_fwhm: float (default 200.0)
                  Full-width half-max of the gaussian core of the line in units of km/s
              lorentzian_fwhm: float (default 500.0)
                  Full-width half-max of the lorentzian wings of the line in units of km/s

background: string (default 'medium') or list-like or numpy.ndarray
  Possible string values are: none, low, medium, and high.  String values trigger the use of
  a notional background model.  If a background spectrum is provided, it is assumed that the
  0th index is the wavelength in microns and the 1st index is the background surface brightness
  in MJy/sr.

calculation: dict
  Set of boolean parameters to toggle the inclusion of different effects and noise parameters in a calculation.
  This section is optional and largely for testing purposes. Do not expect to support this in the UI.

    noise: dict
      Noise components

        darkcurrent: bool
            Dark current
        crs: bool
            Cosmic rays
        rn_correlation: bool
            Correlated read noise
        ffnoise: bool
            Flat-field noise
        readnoise: bool
            Detector read-out noise

    effects: dict
      Effects that can affect the noise or detector response or both

        ipc: bool
            Inter-pixel capacitance
        saturation: bool
            Pixel saturation
        background: bool
            Include background in calculation or not


configuration: dict
  This configuration for the instrument using the following keys:

    instrument: dict
      The instrument configuration parameters

        instrument: string
          for JWST:
            * miri
            * nircam
            * nirspec
            * niriss

          for HST:
            * wfc3 (NOT IMPLEMENTED)

          for WFIRST:
            * wfirstimager
            * wfirstifu

        mode: string
          valid modes:
            * imaging
            * sw_imaging
            * lw_imaging
            * msa
            * mrs
            * soss
            * ifu
            * wfss
            * ssgrism
            * wfgrism
            * lrsslit
            * lrsslitless
            * fixed_slit
            * ami
            * coronagraphy

        filter: string
           (e.g. f070w)

        disperser: string
           (e.g. g235h)

        aperture: string
           (e.g. a200)

        shutter_location: string (only valid for NIRSpec MSA mode)
            Identifier string for slitlet position to use for MSA calculation

        slitlet_shape: list-like  (only valid for NIRSpec MSA mode)
            List of 2-element offsets describing set of shutters to be open. Offsets are from scene center
            in units of shutter spacing.

    detector: dict
      Exposure configuration parameters.

        subarray: string
           full, 64x64, etc.; Instrument-dependent
        readmode: string
           Instrument-dependent
        ngroup: int
           Number of groups
        nint: int
           Number of integrations
        nexp: int
           Number of exposures

    dynamic_scene: boolean
        Toggle whether to allow the size of the scene to expand dynamically to include all configured sources.

    scene_size: float
        Default size of the scene in arcseconds. Used if dynamic_scene is True.

    max_scene_size: float
        Maximum allowable scene_size in arcseconds.

strategy: dict
  Configuration parameters for observing strategy.

    method: string
        Instrument and mode dependent. Currently supported methods are:
            * imagingapphot
            * specapphot
            * coronagraphy
            * ifuapphot
            * ifunodinscene
            * ifunodoffscene
            * msafullapphot
            * soss

        Planned methods that are not yet implemented include:
            imagingoptphot, specoptphot, speclinephot

    units: string  (default: "arcsec")
        Angular units used by the strategy
    target_source: string
        Sent by the UI client, but currently unused by the engine
    target_type: string
        Sent by the UI client, but currently unused by the engine

    The rest of the parameters will be method dependent.  The parameters required
    for **imagingapphot**, **specapphot**, and **ifuapphot** are:

        aperture_size: float
            Size of extraction aperture in "units"
        sky_annulus: list-like of format (float, float)
            The inner and outer radii in "units" of sky region used for background subtraction
        target_xy: list-like of format (float, float)
            X and Y center position of the aperture and sky annulus.

    The parameters required for **ifunodinscene** and **ifunodoffscene** are:

        aperture_size: float
            Size of extraction aperture in "units"
        target_xy: list-like of format (float, float)
            X and Y center position of the aperture and sky annulus.
        dithers: list of dicts with format {'x': <float>, 'y': <float>}
            Dither positions given in "units" from center of the Scene.

    The parameters required for **msafullapphot** are:

        shutter_offset: list-like of format (float, float)
            Offset of shutter pattern from center of scene in "units"
        dithers: list of dicts
            Dither positions and MSA shutter configuration with the following format:
                x: float
                    X position of the central shutter
                y: float
                    Y position of the central shutter
                on_source: list of bool
                    List of booleans denoting whether a shutter should be treated as source or sky.
                    Must be of same length as "slitlet_shape" specified in the instrument configuration.

    The parameters required for **soss** are:
        background_subtraction: boolean
            Toggle whether background subtraction is performed or not.
        order: int
            Specify which order to extract. Can be 1 or 2 with support for 3 forthcoming.

    The parameters required for **coronagraphy** are:

        target_xy: two-element list-like (float, float)
            Position of extraction aperture
        aperture_size: float
            Radius of extraction aperture in 'units'
        sky_annulus: two-element list-like (float, float)
            Inner and outer radii of sky background estimation region
        contrast_azimuth: float
            Azimuth at which to calculate contrast curve
        pointing_error: two-element list-like (float, float)
            Amount to shift occulted source to emulate imperfect pointing
        delta_opd: float
            Change in system OPD
        scene_rotation: float
            Rotation angle to apply to scene
        psf_subtraction_source: Source dict in engine API format
            Definition of source to use for PSF subtraction
        psf_subtraction_xy: two-element list-like (float, float)
            Offset to apply to psf_subtraction_source
        unocculted_xy: two-element list-like (float, float)
            Offset to apply to source to measure contrast between occulted and unocculted observation

    The parameters required for the planned methods will be defined as the methods
    are implemented.

fake_exception: list of strings
    If present, this list is searched for control terms that cause
    perform_calculation to raise exceptions for testing purposes.
    Currently recognized strings are:

        'pyetc':
             raise PyetcException

        'exception':
             raise Exception

    Other strings may be added later to add exceptions or modify
    the details of the exception objects raised.
