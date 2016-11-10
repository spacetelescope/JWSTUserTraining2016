"""
Functions converted from C to python from HEALPix file
chealpix.c

"""
import math

twopi = 2. * math.pi
inv_halfpi = 2. / math.pi
twothird = 2. / 3.


# /*! Returns the remainder of the division \a v1/v2.
#   The result is non-negative.
#    \a v1 can be positive or negative; \a v2 must be positive. */

def fmodulo(v1, v2):
    if v1 >= 0.:
        if (v1 < v2):
            return v1
        else:
            return math.fmod(v1, v2)
    tmp = math.fmod(v1, v2) + v2
    if tmp == v2:
        return 0.
    else:
        return tmp
# /*  return (v1>=0) ? ((v1<v2) ? v1 : fmod(v1,v2)) : (fmod(v1,v2)+v2); */
# 3  }
# *! Returns the remainder of the division \a v1/v2.
#    The result is non-negative.
#    \a v1 can be positive or negative; \a v2 must be positive. */


def imodulo(v1, v2):
    v1 = int(v1)
    v2 = int(v2)
    v = v1 % v2
    if (v >= 0):
        return v
    else:
        return v + v2

# static int isqrt(int v)
#  { return (int)(sqrt(v+0.5)); }


def ang2pix_ring(nside, ra, dec):
    """Converts RA and Dec to HEALPix pixel index.

        Parameters
        ----------
        ra : double
            Right Ascension [radians]
        dec : double
            declination [radians]
        nside : int
            HEALPix scale factor.

        Returns
        -------
        ipix : int
            HEALPix pixel index.
        """

    theta = math.pi / 2. - dec
    phi = ra
    if theta < 0. or theta > math.pi:
        print("theta out of range in ang2pix_ring")
    return ang2pix_ring_z_phi(nside, math.cos(theta), phi)


def ang2pix_ring_z_phi(nside, z, phi):
    za = abs(z)
    tt = fmodulo(phi, twopi) * inv_halfpi  # in [0,4)

    if za <= twothird:  # Equatorial region
        temp1 = nside * (0.5 + tt)
        temp2 = nside * z * 0.75
        jp = int(temp1 - temp2)  # index of  ascending edge line
        jm = int(temp1 + temp2)  # index of descending edge line

        # ring number counted from z=2/3
        ir = nside + 1 + jp - jm  # ; /* in {1,2n+1} */
        kshift = 1 - (ir & 1)  # ; /* kshift=1 if ir even, 0 otherwise */  CHECK THIS &?

        ip = (jp + jm - nside + kshift + 1) / 2  # ; /* in {0,4n-1} */
        ip = imodulo(ip, 4 * nside)

        return nside * (nside - 1) * 2 + (ir - 1) * 4 * nside + ip
    else:  # /* North & South polar caps */
        tp = tt - (int)(tt)  # ;
        tmp = nside * math.sqrt(3 * (1 - za))  # ;

        jp = int(tp * tmp)  # ; /* increasing edge line index */
        jm = int((1.0 - tp) * tmp)  # ; /* decreasing edge line index */

        ir = jp + jm + 1  # ; /* ring number counted from the closest pole */
        ip = int(tt * ir)  # ; /* in {0,4*ir-1} */
        ip = imodulo(ip, 4 * ir)  # ;

        if z > 0.:
            return 2 * ir * (ir - 1) + ip
        else:
            return 12 * nside * nside - 2 * ir * (ir + 1) + ip
