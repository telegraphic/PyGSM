#from .pygsm import GlobalSkyModel

import numpy as np
from scipy.interpolate import interp1d, pchip
import h5py
from astropy import units
import healpy as hp
import pylab as plt
import ephem
from datetime import datetime

from pkg_resources import resource_filename
GSM2016_FILEPATH = resource_filename("pygsm", "gsm2016_components.h5")

kB = 1.38065e-23
C = 2.99792e8
h = 6.62607e-34
T = 2.725
hoverk = h / kB

def K_CMB2MJysr(K_CMB, nu):#in Kelvin and Hz
    B_nu = 2 * (h * nu)* (nu / C)**2 / (np.exp(hoverk * nu / T) - 1)
    conversion_factor = (B_nu * C / nu / T)**2 / 2 * np.exp(hoverk * nu / T) / kB
    return  K_CMB * conversion_factor * 1e20#1e-26 for Jy and 1e6 for MJy

def K_RJ2MJysr(K_RJ, nu):#in Kelvin and Hz
    conversion_factor = 2 * (nu / C)**2 * kB
    return  K_RJ * conversion_factor * 1e20#1e-26 for Jy and 1e6 for MJy


class GlobalSkyModel2016(object):
    """ Global sky model (GSM) class for generating sky models.
    """

    def __init__(self, freq_unit='MHz', unit='MJysr', resolution=0):
        """ Global sky model (GSM) class for generating sky models.

        Upon initialization, the map PCA data are loaded into memory and interpolation
        functions are pre-computed.

        Parameters
        ----------

        Notes
        -----

        """

        if unit not in ['MJysr', 'TCMB', 'TRJ']:
            raise RuntimeError("UNIT ERROR: %s not supported. Only MJysr, TCMB, TRJ are allowed."%unit)
        if resolution == 0:
            resolution = 300

        if resolution < 300:
            nside = 1024
            if freq < 10:
                op_resolution = 48
            else:
                op_resolution = 24
        else:
            nside = 64
            op_resolution = 300

        self.h5 = h5py.File(GSM2016_FILEPATH, "r")
        self.freq_unit = freq_unit

        labels = ['Synchrotron', 'CMB', 'HI', 'Dust1', 'Dust2', 'Free-Free']
        self.n_comp = len(labels)
        self.map_ni_hr = np.array([self.h5['highres_%s_map'%lb][:] for lb in labels])
        self.map_ni_lr = self.h5['lowres_maps']
        self.spec_nf = self.h5['spectra'][:]

        self.pca_map_data = None
        self.interp_comps = None
        #self.update_interpolants()

        self.generated_map_data = None
        self.generated_map_freqs = None

    def generate(self, freqs, resolution=0):
        # convert frequency values into Hz
        freqs = np.array(freqs) * units.Unit(self.freq_unit)
        freqs_ghz = freqs.to('GHz').value

        if isinstance(freqs_ghz, float):
            freqs_ghz = np.array([freqs_ghz])

        try:
            assert np.min(freqs_ghz) >= 0.01
            assert np.max(freqs_ghz) <= 5000
        except AssertionError:
            raise RuntimeError("Frequency values lie outside 10 MHz < f < 5 THz")

        if resolution == 0:
            resolution = 300

        if resolution < 300:
            nside = 1024
            if freq < 10:
                op_resolution = 48
            else:
                op_resolution = 24
        else:
            nside = 64
            op_resolution = 300

        if resolution < 300:
            map_ni = self.map_ni_hr
        else:
            map_ni = self.map_ni_lr

        spec_nf = self.spec_nf
        nfreq = spec_nf.shape[1]

        freq = freqs_ghz[0] * 1e3
        left_index = -1
        for i in range(nfreq - 1):
            if freq >= spec_nf[0, i] and freq <= spec_nf[0, i + 1]:
                left_index = i
                break
        if left_index < 0:
            print "FREQUENCY ERROR: %.2e GHz is outside supported frequency range of %.2e GHz to %.2e GHz."%(freq, spec_nf[0, 0], spec_nf[0, -1])

        interp_spec_nf = np.copy(spec_nf)
        interp_spec_nf[0:2] = np.log10(interp_spec_nf[0:2])
        x1 = interp_spec_nf[0, left_index]
        x2 = interp_spec_nf[0, left_index + 1]
        y1 = interp_spec_nf[1:, left_index]
        y2 = interp_spec_nf[1:, left_index + 1]
        x = np.log10(freq)
        interpolated_vals = (x * (y2 - y1) + x2 * y1 - x1 * y2) / (x2 - x1)
        result = np.sum(10.**interpolated_vals[0] * (interpolated_vals[1:, None] * map_ni), axis=0)

        return result

if __name__ == '__main__':

    g = GlobalSkyModel2016(freq_unit='GHz')

    d = g.generate(1.0)

    import pylab as plt
    import healpy as hp

    hp.mollview(d, nest=True)
    plt.show()

    print g

