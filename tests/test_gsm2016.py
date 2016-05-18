from pygsm.pygsm2016 import GlobalSkyModel2016

g = GlobalSkyModel2016(freq_unit='GHz')

d = g.generate(1.0)

import pylab as plt
import healpy as hp

hp.mollview(d, nest=True)
plt.show()

print g