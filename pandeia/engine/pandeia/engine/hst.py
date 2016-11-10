# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

from .telescope import Telescope
from .instrument import Instrument


class HST(Telescope):

    """
    Currently a dummy class for directory/file discovery, but could eventually contain HST-specific methods
    """
    pass


class HSTInstrument(Instrument):

    """
    Generic HST Instrument class
    """
    def __init__(self, mode=None, config={}, **kwargs):
        telescope = HST()
        Instrument.__init__(self, telescope=telescope, mode=mode, config=config, **kwargs)


class WFC3(HSTInstrument):

    """
    Currently HST WFC3 requires no extra methods beyond those provided by the generic Instrument class
    """
    pass
