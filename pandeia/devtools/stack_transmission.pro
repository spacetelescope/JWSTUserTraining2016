PRO stack_transmission, list, out

spec = MRDFITS(list[0],1)
FOR i=1,N_ELEMENTS(list)-1 DO BEGIN
   spec_new = MRDFITS(list[i],1)
   spec.transmission = spec.transmission * INTERPOL(spec_new.transmission, spec_new.wave, spec.wave) 
ENDFOR

MWRFITS, dummy, out,/CREATE
MWRFITS, {wave:spec.wave,transmission:spec.transmission}, out

END
