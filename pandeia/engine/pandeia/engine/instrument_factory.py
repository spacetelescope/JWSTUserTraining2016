# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

from .custom_exceptions import EngineInputError
from .utils import recursive_subclasses, merge_data
from .instrument import Instrument

# need to import the Instrument subclasses from whereever they're defined
from .jwst import NIRSpec, NIRCam, NIRISS, MIRI
from .wfirst import WFIRSTImager
from .hst import WFC3


def InstrumentFactory(config={}, webapp=False, **kwargs):
    """
    Function to take configuration data and build/return an appropriately configured
    instance of the desired Instrument subclass.

    Parameters
    ----------
    config: dict
        Configuration data in engine API format.  Used to suss out which instrument to configure
        and then passed along.
    **kwargs: keyword/value pairs
        Additional configuration data
    """
    all_config = merge_data(config, dict(**kwargs))
    types = recursive_subclasses(Instrument)
    instruments = [t.__name__.lower() for t in types]
    inst_map = dict(zip(instruments, types))

    # get the instrument name out of the input configuration
    try:
        instrument = all_config['instrument']['instrument']
    except KeyError as e:
        msg = "Must provide Instrument name via engine API compatible dict or instrument=<name> keyword/value pair."
        raise EngineInputError(value=msg)

    # get the instrument mode out of the input configuration
    try:
        mode = all_config['instrument']['mode']
    except KeyError as e:
        msg = "Must provide Instrument mode via engine API compatible dict or mode=<modename> keyword/value pair. (%s)" % repr(e)
        raise EngineInputError(value=msg)

    # make sure requested instrument is one we actually support
    if instrument not in instruments:
        msg = "Instrument %s not supported/implemented." % instrument
        raise EngineInputError(value=msg)
    else:
        cls = inst_map[instrument](mode=mode, config=config, webapp=webapp, **kwargs)
        return cls
