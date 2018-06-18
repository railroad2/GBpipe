#
# This file is part of GBpipe.
#
# GBpipe is a package for GroundBIRD data processing.
#

""" 
Provides observation direction calculation functions 
"""

import healpy as hp
import numpy  as np
import GBparam 
import sys
import time
from astropy.time import Time

R_E2G = hp.Rotator(coord=['C', 'G'])

if sys.version_info < (3,):
    range = xrange

def unixtime2JD(unixtime): # unixtime (float)
    t = Time(unixtime, format='unix') 
    t.format('jd')
    return t.value

def unixtime2LST(unixtime, deg=True):
    t = Time(unixtime, format='unix')
    if (deg==True):
        return t.sidereal_time('apparent', longitude=str(GBparam.lon)+'d').degree
    else: #returns LST in hourangle
        return t.sidereal_time('apparent', longitude=str(GBparam.lon)+'d').value

def JD2LST(jd, deg=True):
    t = Time(jd, format='jd')
    if (deg==True):
        return t.sidereal_time('apparent', longitude=str(GBparam.lon)+'d').degree
    else: #returns LST in hourangle
        return t.sidereal_time('apparent', longitude=str(GBparam.lon)+'d').value

def unixtime2lst(ut):
    import astropy.time
    t = astropy.time.Time(ut, format='unix', location = (GBparam.GBparams.lat, GBparam.GBparam.lon))
    return t.tai.sidereal_time('apparent')

def ut2lst_deg(ut):
    import astropy.time
    t = astropy.time.Time(ut, format='unix', location = (GBparam.GBparam.lat, GBparam.GBparam.lon))
    return t.tai.sidereal_time('apparent').degree

def encoder2ang(enc, south=GBparam.GBparam.encoder_south): # encoder value to angle
    deg = 360.0/8192*(np.array(enc) - south)
    if (not hasattr(deg, '__iter__')):
        deg = [deg]

    for i in range(len(deg)):
        if (deg[i] < 0):
            deg[i] += 360

    return deg

def euler_ZYZ(angles, deg=True): 
    """ Calculates rotation matrix according to the wikipedia convention (extrinsic z-y-z) """
    phi   = np.radians(angles[2]) # alpha
    theta = np.radians(angles[1]) # beta
    psi   = np.radians(angles[0]) # gamma 

    c1 = np.cos(phi)
    c2 = np.cos(theta)
    c3 = np.cos(psi)
    s1 = np.sin(phi)
    s2 = np.sin(theta)
    s3 = np.sin(psi)

    Rtmp = np.array([[ c1*c2*c3-s1*s3,  -c3*s1-c1*c2*s3, c1*s2], 
                     [ c1*s3+c2*c3*s1,   c1*c3-c2*s1*s3, s1*s2], 
                     [-c3*s2,            s2*s3,          c2]])

    if (len(Rtmp.shape)==3):
        R = np.transpose(Rtmp, (2, 0, 1))
    else:
        R = Rtmp

    return R

def xp_coord(angles):
    """ Get x (theta) axis given a direction for parallactic angle calculation """
    theta= angles[0]
    phi  = angles[1]
    cost = np.cos(theta)
    sint = np.sin(theta)
    cosp = np.cos(phi)
    sinp = np.sin(phi)
    xp = np.array([cost * cosp, cost * sinp, -sint])
    #yp = np.array([-sinp, cosp, 0.0])
    return xp #, yp

def parallactic_angle_einsum(ze, deg=True): 
    """ 
    Calculates angle between theta axes given a direction in Equatorial and Galactic coordinates.

    Parameters
    ----------
    ze  : float or float array
        Direction vector in equatorial coordinate.
    deg : bool
        If it is set 'True', the output angle is in degree. 
        Otherwise, the output angle is in radian (default: True).
    
    Returns
    --------
    psi_par : float
        parallactic angle for the given direction
    """
    ze = np.array(ze)
    if (hasattr(ze[0], '__iter__')):
        zg = np.einsum('ij,kj->ki', R_E2G.mat, ze)
        xp = xp_coord(hp.vec2ang(ze)).T
        xe = np.einsum('ij,kj->ki', R_E2G.mat, xp)
        xg = xp_coord(hp.vec2ang(zg)).T
        xcx = np.cross(xe, xg)
        xdx = np.einsum('ij,ij->i', xe, xg)
        idx1 = np.where(xdx > 1)
        xdx[idx1] = 1
        sign = np.sign(np.einsum('ij,ij->i', zg, xcx))
    else:
        zg = R_E2G(ze)
        xp = xp_coord(hp.vec2ang(ze)).flatten()
        xe = R_E2G(xp)
        xg = xp_coord(hp.vec2ang(zg)).flatten()
        xcx = np.cross(xe, xg)
        xdx = np.dot(xe, xg)
        if (xdx > 1): xdx == 1
        sign = np.sign(np.dot(zg, xcx))

    psi_par = np.arccos(xdx) * sign

    return psi_par

def parallactic_angle(ze, deg=True): 
    """ 
    Calculates angle between theta axes given a direction in Equatorial and Galactic coordinates.
    *numpy.einsum()* functions in some products are replaced with numpy.tensordot() for speed. 

    Parameters
    ----------
    ze  : float or float array
        Direction vector in equatorial coordinate.
    deg : bool
        If it is set 'True', the output angle is in degree. 
        Otherwise, the output angle is in radian (default: True).
    
    Returns
    --------
    psi_par : float
        parallactic angle for the given direction
    """
    ze = np.array(ze)
    if (hasattr(ze[0], '__iter__')):
        zg = np.tensordot(R_E2G.mat, ze, axes=(1,1)).T
        xp = xp_coord(hp.vec2ang(ze)).T
        xe = np.tensordot(R_E2G.mat, xp, axes=(1,1)).T
        xg = xp_coord(hp.vec2ang(zg)).T
        xcx = np.cross(xe, xg)
        xdx = np.einsum('ij,ij->i', xe, xg)
        idx1 = np.where(xdx > 1)
        xdx[idx1] = 0
        sign = np.sign(np.einsum('ij,ij->i',zg, xcx))
    else:
        zg = R_E2G(ze)
        xp = xp_coord(hp.vec2ang(ze)).flatten()
        xe = R_E2G(xp)
        xg = xp_coord(hp.vec2ang(zg)).flatten()
        xcx = np.cross(xe, xg)
        xdx = np.dot(xe, xg)
        if xdx > 1: xdx = 1 
        sign = np.sign(np.dot(zg, xcx))

    psi_par = np.arccos(xdx) * sign

    return psi_par

def Rot_matrix_healpix(EL=GBparam.GBparam.EL, AZ=0, LAT=GBparam.GBparam.lat, LST=0, coord='C'): # angles in degree by default
    """Calculates the rotation matrix of GroundBIRD with healpix routine.

    Parameters
    ----------
    EL : 

    Returns
    -------
    """
    
    # rotation matrix calculation (default : Equatorial coordinate 'C')
    r1 = hp.Rotator( (0, 90.-EL, 180.-AZ), eulertype='Y', deg=True) # rotation of GB w.r.t. ground
    if (coord == 'H'):
        rmat = r1.mat
    else: 
        r2 = hp.Rotator( (0, 90.-LAT, LST), eulertype='Y', deg=True) # horizontal to equatorial coordinate. 
        rmat = np.matmul(r2.mat, r1.mat)
        if (coord == 'G'):
            rmat = np.matmul(R_E2G.mat, rmat)

    return rmat

def Rot_matrix(EL=GBparam.GBparam.EL, AZ=0, LAT=GBparam.GBparam.lat, LST=0, coord='C'): # angles in degree by default
    """Computes rotation matrix with euler_ZYZ routine."""
    EL = np.array(EL)
    AZ = np.array(AZ)
    LAT = np.array(LAT)
    LST = np.array(LST)

    r1 = euler_ZYZ((0, 90.-EL, 180.-AZ), deg=True)
    if (coord =='H'):
        rmat = r1
    else:
        r2 = euler_ZYZ((0, 90.-LAT, LST), deg=True)
        rmat = np.matmul(r2, r1)
        if (coord == 'G'):
            r3mat = [R_E2G.mat]
            rmat = np.matmul(r3mat, rmat)

    return rmat

def Rotate(v_arr, rmat=None):
    if (rmat is None):
        rmat = Rot_matrix()

    vp_arr = np.matmul(rmat, np.array(v_arr).T)
    parangle = 0

    #if (len(vp_arr[0])!=1):
    #    vp_arr = vp_arr.T
    #else:
    #    vp_arr = np.ravel(vp_arr.T)

    return vp_arr 

def dir_sky(v_arr):
    rmat = Rot_matrix()
    return Rotate(v_arr, rmat)
    
def dir_pol(v_arr):
    rmat = Rot_matrix()
    return Rotate(v_arr, rmat)

