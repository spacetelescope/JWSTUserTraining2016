Description of ETC API
======================

A simple dict-based API has been implemented as an interface to the JWST ETC engine.
The use of dicts allows for easy serialization of the data going into and out of the engine.
Here is an example of its usage::

  import sys
  sys.path.append("<path_to_jwst_code_base>")
  from perform_calculation import perform_calculation

  output = perform_calculation(input)

where ``input`` is a dict described by :doc:`engine_input_api` and ``output`` is described
by :doc:`engine_output_api`.

There is also a ``fits_report`` function that can be imported from ``perform_calculation``.
It takes the output from ``perform_calculation()`` and converts the 2D and 3D images to
``astropy.io.fits.PrimaryHDU`` format.
