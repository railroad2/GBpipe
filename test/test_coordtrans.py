import healpy as hp
import numpy  as np
import pylab  as plt

def test_compare():
    Nside = 32
    Npix = 12*Nside**2
    m=np.arange(Npix)
    R1 = hp.Rotator(coord=('G', 'C'))
    theta, phi = hp.pix2ang(Nside, m)
    g0, g1 = R1(theta, phi)
    idx1   = hp.ang2pix(Nside, g0, g1)
    m1 = m[idx1]

    h0, h1 = R1.I(theta, phi)
    idx2   = hp.ang2pix(Nside, h0, h1)
    m2 = m[idx2]

    hp.mollview(m,  title='Original', coord='C')
    hp.mollview(m,  title='Transform in Mollview', coord=('C', 'G'))
    hp.mollview(m1, title='Transform by rotation', coord=('G'))
    hp.mollview(m2, title='Transform by inverse rotation', coord=('G'))

    plt.show()
 
def test_compare2():
    Nside = 128
    Npix = 12*Nside**2
    m=np.arange(Npix)
    R1 = hp.Rotator(np.radians((60,30,0)), deg=False)
    theta, phi = hp.pix2ang(Nside, m)
    g0, g1 = R1(theta, phi)
    idx1   = hp.ang2pix(Nside, g0, g1)
    m1 = m[idx1]

    hp.mollview(m, title='Transform in Mollview', rot=(360-30,-60,0))
    hp.mollview(m1, title='Transform by rotation')
    #hp.mollview(m2, title='Transform by inverse rotation')
    plt.show()

def test_compare3(): # compare rotate with matrix and rotator
    Nside = 32 
    Npix 12 * Nside ** 2
    m=np.arange(Npix)
    #R1 = hp.Rotator(np.radians((

def main():
    test_compare()
    #test_compare2()

if __name__=='__main__':
    main()

