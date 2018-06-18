import healpy as hp
import pylab as plt

fname = "m_original.fits"
fname1 = "m_rotation.fits"

m=hp.read_map(fname)
m1=hp.read_map(fname1)

hp.mollview(m)
hp.mollview(m1)
plt.show()
