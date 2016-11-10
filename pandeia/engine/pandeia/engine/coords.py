# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import numpy as np
from photutils.geometry import circular_overlap_grid, elliptical_overlap_grid

from .custom_exceptions import EngineInputError


class Grid(object):

    """
    Generate an x, y grid in a rectangular region, sampled with xsamp and
    ysamp spacings in the x and y directions, respectively.  Origin is at the
    upper-left corner to match numpy and matplotlib convention. astropy.io.fits
    also assumes UL origin by default.

    Parameters
    ----------
    xsamp, ysamp : float, float
        The sampling spacing in the x and y directions.
    nx, ny: int, int
        Number of samples in the x and y directions.
    """

    def __init__(self, xsamp, ysamp, nx, ny):
        self.xsamp = np.abs(xsamp)
        self.ysamp = np.abs(ysamp)
        startx = -(nx - 1) / 2.0 * self.xsamp
        stopx = (nx - 1) / 2.0 * self.xsamp
        starty = -(ny - 1) / 2.0 * self.ysamp
        stopy = (ny - 1) / 2.0 * self.ysamp
        xvals = np.linspace(startx, stopx, num=nx)
        yvals = np.linspace(starty, stopy, num=ny)

        ones = np.ones((ny, nx))
        x = ones * xvals
        y = np.flipud(ones * yvals.reshape(int(ny), 1))

        self.nx = nx
        self.ny = ny
        self.x = x
        self.y = y
        self.row = xvals
        # flip Y axis because we use Y increasing from bottom to top
        self.col = yvals[::-1]

    @property
    def shape(self):
        """
        Provides array-like .shape functionality
        """
        sh = (self.ny, self.nx)
        return sh

    def as_dict(self):
        """
        Create engine API format dict containing the properties of the Grid instance. The half-pixel
        offsets are due to the pixel coordinates being referenced to pixel centers. Thus the edges of the
        Grid are a half-pixel above and to the left of the Grid index origin.
        """
        transform = {}
        transform['x_refpix'] = 0
        transform['x_refval'] = self.x[0][0]
        transform['x_max'] = self.x.max() + 0.5 * self.xsamp
        transform['x_min'] = self.x.min() - 0.5 * self.xsamp
        transform['x_step'] = self.xsamp
        transform['x_size'] = self.x.shape[1]
        transform['y_refpix'] = 0
        transform['y_refval'] = self.y[0][0]
        transform['y_max'] = self.y.max() + 0.5 * self.ysamp
        transform['y_min'] = self.y.min() - 0.5 * self.ysamp
        transform['y_step'] = -self.ysamp
        transform['y_size'] = self.y.shape[0]
        return transform

    def wcs_info(self):
        """
        Define coordinate transform in WCS header format.

        Returns
        -------
        header: dict
            WCS keys defining coordinate transform for the 2 spatial axes
        """
        t = self.as_dict()
        header = {
            'ctype1': 'X offset',
            'crpix1': 1,
            'crval1': t['x_min'],
            'cdelt1': t['x_step'],
            'cunit1': 'arcsec',
            'cname1': 'X',
            'ctype2': 'Y offset',
            'crpix2': 1,
            'crval2': t['y_min'],
            'cdelt2': -t['y_step'],
            'cunit2': 'arcsec',
            'cname2': 'Y',
        }
        return header

    def __len__(self):
        """
        Provide array-like len() functionality
        """
        l = self.x.shape[0] * self.x.shape[1]
        return l

    def __getitem__(self, val):
        """
        Provide array-like indexing functionality.

        Parameters
        ----------
        val : slice object
            Valid python index or slice() specification.

        Returns
        -------
        pos : [y, x] list of floats or arrays
            Elements are floats in the case of specific indicies and arrays in
            the cases of slice()'s.
        """
        section = [self.y[val], self.x[val]]
        return section

    def bounds(self):
        return {'xmin': np.min(self.x),
                'xmax': np.max(self.x),
                'ymin': np.min(self.y),
                'ymax': np.max(self.y)}

    def dist(self, xcen=0.0, ycen=0.0):
        """
        Return a distance array where each element contains its distance from
        the center of the grid.
        """
        d = np.sqrt((self.x - xcen) ** 2 + (self.y - ycen) ** 2)
        return d

    def world_to_image(self, y, x):
        """
        Return fractional index coordinates yi,xi corresponding to given y,x position

        Parameters
        ----------
        y, x : float
            y, x position in world coordinates

        Returns
        -------
        pos : [yi, xi] list of floats
            yi, xi image coordinates that correspond to given position
        """
        t = self.as_dict()
        if x < t['x_min'] or x > t['x_max']:
            raise ValueError("X outside of Grid bounds.")
        if y < t['y_min'] or y > t['y_max']:
            raise ValueError("Y outside of Grid bounds.")
        x_image = x / self.xsamp + (self.shape[1] - 1) / 2.0
        # there's a minus sign here because y indices increase from up to down,
        # but y coordinates increase down to up.
        y_image = -y / self.ysamp + (self.shape[0] - 1) / 2.0
        pos = [y_image, x_image]
        return pos

    def world_to_index(self, y, x):
        """
        Return nearest index pair to given y,x position

        Parameters
        ----------
        y, x : float
            y,x position in world coordinates

        Returns
        -------
        pos : [y, x] list of ints
            y,x indices that most closely correspond to given position
        """
        y_index, x_index = self.world_to_image(y, x)
        y_index, x_index = self._index_bound(int(round(y_index)), self.ny), self._index_bound(int(round(x_index)), self.nx)
        pos = [y_index, x_index]
        return pos

    def _index_bound(self, index, nindex):
        """
        Apply bounds to indicies so the nearest index is either the first or last index for cases
        at edge or beyond Grid.

        Parameters
        ----------
        index: int
            Index to bounds-check
        nindex: int
            Length of axis being indexed

        Returns
        -------
        index: int
            Index with boundaries applied
        """
        if index < 0:
            index = 0
        if index > nindex - 1:
            index = nindex - 1
        return index

    def shift_rotate(self, yoff, xoff, rot):
        """
        Return shifted/rotated (y, x) given offsets (yoff, xoff) and rotation, rot (degrees)

        Parameters
        ----------
        yoff, xoff: float
            yoff, xoff offsets in world coordinates
        rot: float
            rotation angle in degrees

        Returns
        -------
        ysh_rot, xsh_rot: 2D numpy arrays
            rotated and shifted copies of Grid.x and Grid.y
        """
        pa_radians = np.pi * rot / 180.0
        xsh = self.x - xoff
        ysh = self.y - yoff
        xsh_rot = xsh * np.cos(pa_radians) + ysh * np.sin(pa_radians)
        ysh_rot = -xsh * np.sin(pa_radians) + ysh * np.cos(pa_radians)
        return ysh_rot, xsh_rot

    def get_aperture(self):
        """
        Return Grid parameters as an aperture specification
        """
        ap = {'width': self.nx * self.xsamp, 'height': self.ny * self.ysamp, 'offset': (0, 0)}
        return ap

    def rectangular_mask(self, width=None, height=None, xoff=0.0, yoff=0.0, transparency=1.0):
        """
        Define a rectangular mask within the Grid of specfied width and height and offset by xoff/yoff.
        Pixels partially covered by mask are properly weighted by the area subtended by the mask.

        Parameters
        ----------
        width: float
            Width of the mask. If None, then use the full width of the Grid.
        height: float
            Height of the mask. If None, then use the full height of the Grid.
        xoff: float
            X offset of the mask
        yoff: float
            Y offset of the mask
        transparency: float (default: 1.0)
            Transparency of the mask

        Returns
        -------
        mask: 2D np.ndarray
            2D mask image
        """
        if transparency < 0.0 or transparency > 1.0:
            msg = "Mask transparency, %f, must be in the range of 0.0 (fully opaque) to 1.0 (fully clear)." % transparency
            raise EngineInputError(value=msg)

        minx = self.row.min()
        maxx = self.row.max()
        xcen = (minx + maxx) / 2.0
        miny = self.col.min()
        maxy = self.col.max()
        ycen = (miny + maxy) / 2.0

        # if width and/or height are not specified, set them encompass the size of the Grid.
        # add two pixels to make sure the masking does not interpolate the edges.
        if width is None:
            width = maxx - minx + 2 * self.xsamp
        if height is None:
            height = maxy - miny + 2 * self.xsamp

        mask_x = np.array(
            [
                minx,
                xcen - (width + self.xsamp) / 2.0 + xoff,
                xcen - (width - self.xsamp) / 2.0 + xoff,
                xcen + (width - self.xsamp) / 2.0 + xoff,
                xcen + (width + self.xsamp) / 2.0 + xoff,
                maxx
            ]
        )
        mask_y = np.array(
            [
                miny,
                ycen - (height + self.ysamp) / 2.0 + yoff,
                ycen - (height - self.ysamp) / 2.0 + yoff,
                ycen + (height - self.ysamp) / 2.0 + yoff,
                ycen + (height + self.ysamp) / 2.0 + yoff,
                maxy
            ]
        )

        mask_pix = np.array([0.0, 0.0, 1.0, 1.0, 0.0, 0.0])
        x_mask_int = np.interp(self.row, mask_x, mask_pix)
        y_mask_int = np.interp(self.col, mask_y, mask_pix)
        x_mask = x_mask_int.reshape((1, self.nx))
        y_mask = y_mask_int.reshape((self.ny, 1))
        mask = y_mask * x_mask * transparency
        return mask

    def circular_mask(self, radius, xoff=0.0, yoff=0.0, use_exact=1, subsampling=1, transparency=1.0):
        """
        Define a circular mask within the grid of a specified radius and offset by xoff/yoff.
        Use photutils.geometry to create masks that calculate the fraction of a pixel subtended
        by the circular mask.  use_exact=1 performs the exact geometric calculation while
        use_exact=0 will sub-sample the pixels by the specfied subsampling parameter
        to estimate the subtended area.

        Parameters
        ----------
        radius: float
            Radius of the circular mask
        xoff: float
            X position of the center of the mask
        yoff: float
            Y position of the center of the mask
        use_exact: int (default: 1)
            If 1, then use exact geometrical calculation to weight pixels partially covered by mask.
            This parameter is passed directly on to photutils.geometry.circular_overlap_grid()
        subsampling: int (default: 1)
            If use_exact=0, then subsample pixels by this factor to calculate area of each pixel subtended
            by the mask. The default value of 1 will force no partial pixels to be included in the mask.
            This parameter is passed directly on to photutils.geometry.circular_overlap_grid()
        transparency: float (default: 1.0)
            Transparency of the mask
        Returns
        -------
        mask: 2D np.ndarray
            2D mask image
        """
        if transparency < 0.0 or transparency > 1.0:
            msg = "Mask transparency, %f, must be in the range of 0.0 (fully opaque) to 1.0 (fully clear)." % transparency
            raise EngineInputError(value=msg)

        t = self.as_dict()
        # we need to use flipud because we use an origin in the UL corner of an image
        # while photutils uses the LL corner.
        mask = np.flipud(
            circular_overlap_grid(
                t['x_min'] - xoff,
                t['x_max'] - xoff,
                t['y_min'] - yoff,
                t['y_max'] - yoff,
                t['x_size'],
                t['y_size'],
                radius,
                use_exact,
                subsampling
            )
        )
        mask *= transparency
        return mask

    def elliptical_mask(self, major, minor, pa=0.0, xoff=0.0, yoff=0.0, use_exact=1, subsampling=1, transparency=1.0):
        """
        Define an elliptical mask within the grid with specified major/minor axes, position angle,
        and offset from center by xoff/yoff. Use photutils.geometry to create masks that calculate
        the fraction of a pixel subtended by the mask.  use_exact=1 performs the exact geometric
        calculation while use_exact=0 will sub-sample the pixels by the specfied subsampling parameter
        to estimate the subtended area.

        Parameters
        ----------
        major: float
            Semi-major axis of the elliptical mask
        minor: float
            Semi-minor axis of the mask
        pa: float
            Position angle of the ellipse in degrees measured clockwise (positive in +X direction)
        xoff: float
            X position of the center of the mask
        yoff: float
            Y position of the center of the mask
        use_exact: int (default: 1)
            If 1, then use exact geometrical calculation to weight pixels partially covered by mask.
            This parameter is passed directly on to photutils.geometry.circular_overlap_grid()
        subsampling: int (default: 1)
            If use_exact=0, then subsample pixels by this factor to calculate area of each pixel subtended
            by the mask. The default value of 1 will force no partial pixels to be included in the mask.
            This parameter is passed directly on to photutils.geometry.circular_overlap_grid()
        transparency: float (default: 1.0)
            Transparency of the mask

        Returns
        -------
        mask: 2D np.ndarray
            2D mask image
        """
        if transparency < 0.0 or transparency > 1.0:
            msg = "Mask transparency, %f, must be in the range of 0.0 (fully opaque) to 1.0 (fully clear)." % transparency
            raise EngineInputError(value=msg)

        theta = np.pi * pa / 180.0  # photutils uses angles in radians

        t = self.as_dict()
        # we need to use flipud because we use an origin in the UL corner of an image
        # while photutils uses the LL corner.
        mask = np.flipud(
            elliptical_overlap_grid(
                t['x_min'] - xoff,
                t['x_max'] - xoff,
                t['y_min'] - yoff,
                t['y_max'] - yoff,
                t['x_size'],
                t['y_size'],
                major,
                minor,
                theta,
                use_exact,
                subsampling
            )
        )
        mask *= transparency
        return mask

    def point_source(self, xoff=0.0, yoff=0.0):
        """
        Define a point source within the grid at position xoff, yoff.  Sub-pixel positioning
        of the point source is performed by weighting the 4 nearest pixels.

        Parameters
        ----------
        xoff: float
            X position of point source
        yoff: float
            Y position of point source

        Returns
        -------
        im: 2D np.ndarray
            Image containing point source with flux scaled to 1.0.
        """
        im = np.zeros(self.shape)

        # normally it's an error to specify a position off the Grid, but in this case it just
        # means there's no flux being added onto the image. in cases such as dithering, putting
        # the point source off the Grid is intentional and shouldn't generate errors or warnings.
        try:
            yi, xi = self.world_to_image(yoff, xoff)
        except ValueError as e:
            return im

        rem_x, int_x = np.modf(xi)
        rem_y, int_y = np.modf(yi)

        # check for image edges and handle weights correctly if point source straddles an edge
        if rem_x < 0.0:
            x_shift = -1
            rem_x = 1.0 + rem_x
        else:
            x_shift = 0

        if rem_y < 0.0:
            y_shift = -1
            rem_y = 1.0 + rem_y
        else:
            y_shift = 0

        # here we handle the cases where a point source is at the edge of the Grid so that some of the
        # flux is in pixels on the Grid, while the rest is off the Grid.
        ly = int(int_y + y_shift)
        uy = int(ly + 1)
        lx = int(int_x + x_shift)
        ux = int(lx + 1)

        # if a pixel is None, then it's off the Grid and we don't try to set it.
        if ly < 0:
            ly = None
        if uy == self.ny:
            uy = None
        if lx < 0:
            lx = None
        if ux == self.nx:
            ux = None

        if ly is not None and lx is not None:
            im[ly, lx] += (1.0 - rem_x) * (1.0 - rem_y)
        if uy is not None and lx is not None:
            im[uy, lx] += (1.0 - rem_x) * rem_y
        if ly is not None and ux is not None:
            im[ly, ux] += rem_x * (1.0 - rem_y)
        if uy is not None and ux is not None:
            im[uy, ux] += rem_x * rem_y

        return im


class IrregularGrid(Grid):

    """
    Create Grid from pre-specified rows and cols.
    """

    def __init__(self, yvals, xvals):
        ny = yvals.size
        nx = xvals.size
        ones = np.ones((ny, nx))
        y = ones * yvals.reshape(ny, 1)
        x = ones * xvals

        self.nx = nx
        self.ny = ny
        self.x = x
        self.y = y
        self.row = xvals
        self.col = yvals
        if len(xvals) > 1:
            self.xsamp = np.abs(xvals[1] - xvals[0])
        if len(yvals) > 1:
            self.ysamp = np.abs(yvals[1] - yvals[0])
