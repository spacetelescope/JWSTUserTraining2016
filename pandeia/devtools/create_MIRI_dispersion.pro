;Reference FM MRS-OPT-08 Wavelength characterization, Martinez-Galarza
;& Meyer 2011 and FM MRS-OPT-07 spectral resolution, Martinez-Galarza
;2011.
;

FUNCTION disp, z, coeff

  RETURN, coeff[0]*z^2.+coeff[1]*z+coeff[2]

END

PRO create_MIRI_dispersion
  PATH = '/Users/pontoppi/STSCI/SoftwareDevelopment/scamp/refdata/MIRI/dispersion/'

  z = findgen(1500) ; Significantly more than 1024, presumably because the spectra are curved?
  
  sampling = [1.9,2.0,2.0,2.2,2.2,2.2,2.4,3.1,3.0,3.3,3.4,4.1]

  coeff_1A = [1e-11,7.67e-4,4.872]
  coeff_1B = [1e-9, 8.62e-4, 5.571]
  coeff_1C = [1e-9, 9.99e-4, 6.424]
  
  coeff_2A = [1e-9, 11.1e-4, 7.270]
  coeff_2B = [1e-9, 13.7e-4, 8.691]
  coeff_2C = [1e-10, 15.4e-4, 9.869]
  
  coeff_3A = [10.^(-9.5), 17.7e-4, 11.405]
  coeff_3B = [1e-10, 20.3e-4, 13.125]
  coeff_3C = [1e-9, 23.8e-4, 15.305]
  
  coeff_4A = [1e-8, 29.9e-4, 17.71]
  coeff_4B = [1e-9, 29.6e-4, 20.495]
  coeff_4C = [1e-8, 28.75e-4, 23.705]
  
  coeffs = [[coeff_1A], [coeff_1B], [coeff_1C],$
            [coeff_2A], [coeff_2B], [coeff_2C],$
            [coeff_3A], [coeff_3B], [coeff_3C],$
            [coeff_4A], [coeff_4B], [coeff_4C]]
  
  filenames = ['ch1_short_disp.fits', 'ch1_medium_disp.fits', 'ch1_long_disp.fits',$
               'ch2_short_disp.fits', 'ch2_medium_disp.fits', 'ch2_long_disp.fits',$
               'ch3_short_disp.fits', 'ch3_medium_disp.fits', 'ch3_long_disp.fits',$
               'ch4_short_disp.fits', 'ch4_medium_disp.fits', 'ch4_long_disp.fits']
  

  FOR i=0,11 DO BEGIN
     disp = disp(z,coeffs[*,i])
     dlds = deriv(z,disp)
     R    = disp/dlds/sampling[i]
     MWRFITS, dummy, PATH+filenames[i], /CREATE
     MWRFITS, {pixels:z,wave:disp,dlds:dlds,R:R, wave_unit:'microns'}, PATH+filenames[i]
  END

END

