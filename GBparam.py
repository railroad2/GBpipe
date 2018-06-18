from __future__ import print_function
import numpy as np
import os.path
import sys
if (sys.version_info > (3,)):
    import configparser as cparser
else:
    import ConfigParser as cparser
from astropy.utils.console import color_print

def perror(message, *args, **kwargs):
    color_print(' ERROR: '+message, 'red', *args, file=sys.stderr, **kwargs)

def pwarning(message, *args, **kwargs):
    color_print(' WARNING: '+message, 'yellow', *args, file=sys.stderr, **kwargs)

class GBparam:
    """ Stores parameters of GroundBIRD.
    All angles are in degree.

    Members
    -------
    tilt : float
        Tilt angle of the telescope
        default = 30
    EL : float
        Elevation of apperture center
        90 - tilt (in degree)
        default = 60
    lat : float
        Lattitude of the telescope
        default = lat_canary =  28d16m7s N =  28.268611
    lon : float
        Longitude of the telescope
        default = lon_canary = 16d36m20s W = -16.605555
    omega_gb : float
        Angular speed of the telescope
        default = 120 degree/s (20 rpm)
    omega_earth : float
        Angular speed of the Earth' rotation
        1 rotation in 1 sidereal day (~ 86164s)
        default = 360 / 86164
    encoder_south : integer
        Encoder value of the South
    fname_pixel : string
        Name of the pixel information file.
        If it is not defined or the file does not exist, 
        dummy data with one pixel is used. 
    """
    tilt          = 30.
    EL            = 90. - 30.
    lat           = (28. + 16./60 +  7./3600)
    lon           = -1 * (16. + 36./60 + 20./3600)
    fsample       = 1200         # sample/s
    omega_gb      = 360. / 3     # degree/s
    omega_earth   = 360. / 86164 # degree/s
    encoder_south = 3185         # at KCH
    fname_pixel   = 'pixelinfo.dat' 
    pixinfo       = np.array([(0, 0, 0, 0., 0., 0., 0., 0., 0.)], 
                            dtype=[('Npix','int32'), ('Nmodule', 'int32'), ('Npix_mod', 'int32'), 
                                   ('X_fc', 'float64'), ('Y_fc', 'float64'), 
                                   ('theta', 'float64'), ('phi', 'float64'), 
                                   ('psi_fc', 'float64'), ('psi_far', 'float64')])

    def __init__(self, fname='default.ini'):
        if (not os.path.isfile(fname)):
            pwarning('File "%s" does not exist. Using internal values.' % fname)
        else:
            self.load_settings(fname) 

        self.load_pixelInfo(self.fname_pixel)

    def get_option(self, cp, sect, opt):
        try:
            return cp.get(sect, opt)
        except cparser.NoOptionError:
            pwarning('No option {0} in {1}'.format(opt, sect))
            return getattr(self, opt)
                 
    def load_settings(self, fname):
        cp = cparser.RawConfigParser()
        cpath   = fname
        cp.read(cpath)

        self.tilt           = float(self.get_option(cp, 'GB', 'tilt'))
        self.lat            = float(self.get_option(cp, 'GB', 'lat'))
        self.lon            = float(self.get_option(cp, 'GB', 'lon'))
        self.rot_speed      = float(self.get_option(cp, 'GB', 'rot_speed'))
        self.fsample        = int(self.get_option(cp, 'GB', 'fsample'))
        self.encoder_south  = int(self.get_option(cp, 'GB', 'encoder_south'))
        self.fname_pixel    = self.get_option(cp, 'GB', 'fname_pixel')
        self.omega_earth    = float(eval(self.get_option(cp, 'others', 'omega_earth')))

        self.EL             = 90. - self.tilt
        self.omega_gb       = 360. * self.rot_speed/60.

    def load_pixelInfo(self, fname=fname_pixel):
        try:
            if (fname[-3:-1]=='csv'):
                self.pixinfo = np.genfromtxt(fname, names=True, delimiter=',')
            else:
                self.pixinfo = np.genfromtxt(fname, names=True)
        except IOError:
            pwarning('Pixel information file "%s" does not exist. Using dummy focalplane.' % fname)

    def show_parameters(self):
        print('-'*50)
        print('GroundBIRD parameters')
        print('-'*50)
        print('Tilt      = ', self.tilt, '(deg)')
        print('Elevation = ', self.EL, '(deg)')
        print('Latitude  = ', self.lat, '(deg)')
        print('Longitude = ', self.lon, '(deg)')
        print('sampling frequency = ', self.fsample, '(sample/s)')
        print('rotation speed = ', self.omega_gb, '(deg/s)')
        print('earth rotation speed = ', self.omega_earth, '(deg/s)')
        print('encoder of South direction = ', self.encoder_south)
        print('-'*50)
        print()

    def show_pixelInfo(self):
        print('-'*50)
        print('GroundBIRD pixel information')
        print('-'*50)
        print(self.pixinfo.dtype.names)
        print(self.pixinfo)
        print('-'*50)
        print()
    
def test_parser():
    par = GBparameters()
    par.show_pixelInfo()
    par.show_parameters()

if __name__=='__main__':
    test_parser()
