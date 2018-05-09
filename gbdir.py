# 
# This file is part of GBpipe.
#
# GBpipe is a package for GroundBIRD data processing.
#

""" Provides observation direction calculation functions 
"""

import healpy as hp
import numpy  as np

def calcJD():
    pass 

def calcLST():
    pass

class GBangles(): # angles in degree
    """ Stores the angle properties of GroundBIRD.
    All angles are in degree.

    Members
    -------
    tilt : float
        tilt angle of the telescope
        default = 30
    zen_EL : float
        Elevation of zenith 
        90 - tilt (in degree)
        default = 60
    lat : float
        lattitude of the telescope
        default = 90 - lat_canary = 61.73138888 where lat_canary = 28d16m7s
    lon : float
        longitude of the telescope
        default = lon_canary = -16.605555 where lon_canary = 16d36m20s
    omega_gb : float
        angular speed of the telescope
        default = 120 degree/s (20 rpm)
    omega_earth : float
        angular speed of the Earth' rotation
        1 rotation in 1 sidereal day (~ 86164s)
        default = 360 / 86164 
    """
    tilt = 30.
    zen_EL   = 90.-30.
    lat  = 90 - (28. + 16./60 +  7./3600)
    lon  = -1 * (16. + 36./60 + 20./3600)
    omega_gb    = 360. / 3
    omega_earth = 360. / 86164

def GBrotmat(EL=GBangles.zen_EL, AZ=0, LAT=GBangles.lat, LST=0, coord='C'): # angles in degree
    """Calculates the rotation matrix of GroundBIRD. 

    Parameters
    ----------
    EL

    Returns
    -------
    """
    r1 = hp.Rotator( (0, 90.-EL, 180.-AZ), eulertype='Y', deg=True)   # rotation of GB w.r.t. ground
    r2 = hp.Rotator( (0, LAT, LST), eulertype='Y', deg=True) # horizontal to equatorial coordinate.
    res = np.matmul(r2.mat, r1.mat)

    if (coord == 'G'):
        r3 = hp.Rotator( coord=['C', 'G'] ) # equatorial to galactic coordinate. 
        res = np.matmul(r3.mat, res)

    return res


#modules for tests
def test_GBrotmat():
    print(GBrotmat())

if __name__=='__main__':
    test_GBrotmat()
