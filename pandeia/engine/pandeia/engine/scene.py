# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import division, absolute_import

import numpy as np

from . import source as src
from .custom_exceptions import EngineInputError


class Scene(object):

    """
    Container class describing a scene as composed of a list of Sources.

    Parameters
    ----------
    input : dict or list (optional, default None)
        Dictionary input must contain a 'Scene' section concordant with the engine input API.
        Can also provide a pre-made list of source.Source objects. If no input supplied,
        then default to single default source.Source().

    Attributes
    ----------
    sources : list
        The list of sources that compose the scene
    offset : list-like (optional)
        Offset to apply to sources upon instantiation given as [x, y]

    Methods
    -------
    offset(x=<float, default 0.0>, y=<float, default 0.0>) :
        Apply x and/or y offset to all sources in scene
    get_size():
        Get the full extent of the scene.
    """

    def __init__(self, input=None, offset=None, webapp=False, **kwargs):
        self.webapp = webapp
        self.warnings = {}
        self.sources = []
        if isinstance(input, dict):
            try:
                # if we get a full engine API dict, fish the scene out
                if 'scene' in input['scene']:
                    self._from_list(input['scene'])
                else:
                    # if we're a dict, but not a full engine API one, try adding as a source
                    self.add_sources(input)
            except KeyError as e:
                message = "Input dict must have source information in engine API format: %s" % e
                raise EngineInputError(value=message)
        elif isinstance(input, list):
            self._from_list(input)
        else:
            self.sources = [src.Source()]

        if offset is not None:
            if isinstance(offset, list) or isinstance(offset, tuple):
                if len(offset) == 2:
                    self.offset({'x': offset[0], 'y': offset[1]})
                elif isinstance(offset, dict):
                    self.offset(offset)
                else:
                    message = "Offset needs to be list-like of format [x, y] "
                    message += "or a dict of format {'x': x, 'y': y}"
                    raise EngineInputError(value=message)
            else:
                message = "If specified, offset needs to be list-like of format [x, y] "
                message += "or a dict of format {'x': x, 'y': y}."
                raise EngineInputError(value=message)

        for s in self.sources:
            self.warnings.update(s.warnings)

    def _from_list(self, source_list):
        """
        Add sources to a scene from a list. Can consist of src.Source() instances or dicts that
        can be used to instantiate a src.Source() or a mixture of the two.

        Parameters
        ----------
        source_list: list
            List consisting of src.Source() instances or dicts for configuring a src.Source()
        """
        for s in source_list:
            if isinstance(s, src.Source):
                self.sources.append(s)
            if isinstance(s, dict):
                self.sources.append(src.Source(config=s, webapp=self.webapp))

    def add_sources(self, s):
        """
        Add a source to a scene. Can be either src.Source instance or a dict containing
        configuration information for a src.Source.

        Parameters
        ----------
        s: src.Source or dict
            A source instance or configuration information for creating source
        """
        if isinstance(s, src.Source):
            self.sources.append(s)
        elif isinstance(s, dict):
            self.sources.append(src.Source(config=s, webapp=self.webapp))
        elif isinstance(s, list):
            self._from_list(s)
        else:
            message = "Must provide either Source instance(s) or a configuration dict(s) to create them."
            raise EngineInputError(value=message)

    def offset(self, dither):
        """
        Apply an offset to all sources within the Scene instance.

        Parameters
        ----------
        dither : dict - {'x': x, 'y': y}
            x : float (default 0.0)
                offset in the X direction
            y : float (default 0.0)
                offset in the Y direction
        """
        for i, s in enumerate(self.sources):
            self.sources[i].position.update(
                {
                    'x_offset': s.position['x_offset'] + dither['x'],
                    'y_offset': s.position['y_offset'] + dither['y']
                }
            )

    def rotate(self, angle=0.0):
        """
        Apply a rotation to the Scene

        Parameters
        ----------
        angle: float
            Rotation angle to apply to Scene. Positive angle rotates Scene CCW.
        """
        angle = np.radians(angle)
        for i, s in enumerate(self.sources):
            x = s.position['x_offset']
            y = s.position['y_offset']
            self.sources[i].position.update(
                {
                    'x_offset': x * np.cos(angle) - y * np.sin(angle),
                    'y_offset': x * np.sin(angle) + y * np.cos(angle)
                }
            )

    def set_pixsamp(self, xsamp=1.0, ysamp=1.0):
        """
        Set the pixel sampling to apply to the scene

        Parameters
        ----------
        xsamp: float
            Pixel size along the X axis
        ysamp: float
            Pixel size along the Y axis
        """
        for s in self.sources:
            s.shape['xsamp'] = xsamp
            s.shape['ysamp'] = ysamp

    def get_size(self):
        """
        Get the relative extent of the scene based on the location of the sources in it.

        Parameters
        ----------
        None

        Returns
        -------
        size: float
            Size of the scene
        """

        # The scene size is the maximum (or minimum) x or y offset.
        # It is multiplied by 2 because the 0 position is referenced to the center of the FOV.
        size = 2.0 * np.max([np.max(np.abs([s.position['x_offset'], s.position['y_offset']])) for s in self.sources])

        return size

    def get_min_size(self):
        """
        Get the minimum scene size required to encompass at least one of the sources.

        Parameters
        ----------
        None

        Returns
        -------
        min_size: float
            Minimum size of the scene
        """

        # The scene size is the maximum (or minimum) x or y offset.
        # It is multiplied by 2 because the 0 position is referenced to the center of the FOV.
        min_size = 2.0 * np.min([np.max(np.abs([s.position['x_offset'], s.position['y_offset']])) for s in self.sources])

        return min_size
