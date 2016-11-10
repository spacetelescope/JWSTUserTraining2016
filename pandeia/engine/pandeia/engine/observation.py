# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

from .config import DefaultConfig
from .custom_exceptions import EngineInputError
from .scene import Scene
from .strategy import StrategyFactory


class Observation(DefaultConfig):

    """
    This class collects the pieces required to make an ETC simulated observation.
    """

    def __init__(self, scene=None, instrument=None, strategy=None, config={}, webapp=False, **kwargs):
        DefaultConfig.__init__(self, config=config, **kwargs)

        # Initialize the observation with the default source and the default strategy
        # relevant for the instrument/mode combo. The user can also pass an initializing source.
        if scene is None:
            self.scene = Scene(webapp=webapp)
        else:
            self.scene = scene

        # make sure we're given a properly configured Instrument
        if instrument is not None:
            self.instrument = instrument
            self.mode = self.instrument.mode
        else:
            msg = "Must provide configured Instrument subclass."
            raise EngineInputError(value=msg)

        # use strategy we're given or fall back to default for the given instrument
        if strategy is not None:
            self.strategy = strategy
        else:
            self.strategy = StrategyFactory(self.instrument, webapp=webapp)

    def get_random_seed(self):
        """
        Return the configured seed for the random number generator.
        """
        return self.random_seed
