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


def rotate_map(hmap, rot_theta, rot_phi, nest=True):
    nside = hp.npix2nside(len(hmap))

    # Get theta, phi for non-rotated map
    t, p = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), nest= nest)  # theta, phi

    # Define a rotator
    r = hp.Rotator(deg=False, rot=[rot_phi, rot_theta])

    # Get theta, phi under rotated co-ordinates
    trot, prot = r(t, p)

    # Inerpolate map onto these co-ordinates
    rot_map = hp.get_interp_val(hmap, trot, prot, nest= nest)

    return rot_map


class GlobalSkyModel2016(object):
    """ Global sky model (GSM) class for generating sky models.
    """

    def __init__(self, freq_unit='MHz', unit='TCMB', resolution='hi', theta_rot=0, phi_rot=0):
        """ Global sky model (GSM) class for generating sky models.

        Upon initialization, the map PCA data are loaded into memory and interpolation
        functions are pre-computed.

        Parameters
        ----------
        freq_unit: 'Hz', 'MHz', or 'GHz'
            Unit of frequency. Defaults to 'MHz'.
        unit: 'MJysr', 'TCMB', 'TRJ'
            Unit of output data. MJy/Steradian, T_CMB in Kelvin, or T_RJ.
        resolution: 'hi' or 'low'
            Resolution of output map. Either 300 arcmin (low) or 24 arcmin (hi).
            For frequencies under 10 GHz, output is 48 arcmin.

        Notes
        -----

        """

        if unit not in ['MJysr', 'TCMB', 'TRJ']:
            raise RuntimeError("UNIT ERROR: %s not supported. Only MJysr, TCMB, TRJ are allowed." % unit)

        if resolution.lower() in ('hi', 'high', 'h'):
            resolution = 'hi'
        elif resolution.lower() in ('low', 'lo', 'l'):
            resolution = 'low'
        else:
            raise RuntimeError("RESOLUTION ERROR: Must be either hi or low, not %s" % resolution)

        self.h5 = h5py.File(GSM2016_FILEPATH, "r")
        self.freq_unit = freq_unit
        self.unit = unit
        self.resolution = resolution

        # Map data to load
        labels = ['Synchrotron', 'CMB', 'HI', 'Dust1', 'Dust2', 'Free-Free']
        self.n_comp = len(labels)

        if resolution=='hi':
            self.map_ni = np.array([self.h5['highres_%s_map'%lb][:] for lb in labels])
        else:
            self.map_ni = np.array(self.h5['lowres_maps'])

        self.spec_nf = self.h5['spectra'][:]

        if theta_rot or phi_rot:
            for i,map in enumerate(self.map_ni):
                self.map_ni[i] = rotate_map(map, theta_rot, phi_rot, nest=True)

        self.generated_map_data = None
        self.generated_map_freqs = None

    def generate(self, freqs):
        """ Generate a global sky model at a given frequency or frequencies

        Parameters
        ----------
        freqs: float or np.array
            Frequency for which to return GSM model

        Returns
        -------
        gsm: np.array
            Global sky model in healpix format, with NSIDE=1024. Output map
            is in galactic coordinates, ring format.

        """

        # convert frequency values into Hz
        freqs = np.array(freqs) * units.Unit(self.freq_unit)
        freqs_ghz = freqs.to('GHz').value

        if isinstance(freqs_ghz, float):
            freqs_ghz = np.array([freqs_ghz])

        try:
            assert np.min(freqs_ghz) >= 0.01
            assert np.max(freqs_ghz) <= 5000
        except AssertionError:
            raise RuntimeError("Frequency values lie outside 10 MHz < f < 5 THz: %s")

        map_ni = self.map_ni
        # if self.resolution == 'hi':
        #     map_ni = self.map_ni_hr
        # else:
        #     map_ni = self.map_ni_lr

        spec_nf = self.spec_nf
        nfreq = spec_nf.shape[1]

        output = np.zeros((len(freqs_ghz), map_ni.shape[1]), dtype='float32')
        for ifreq, freq in enumerate(freqs_ghz):

            left_index = -1
            for i in range(nfreq - 1):
                if spec_nf[0, i] <= freq <= spec_nf[0, i + 1]:
                    left_index = i
                    break

            # Do the interpolation
            interp_spec_nf = np.copy(spec_nf)
            interp_spec_nf[0:2] = np.log10(interp_spec_nf[0:2])
            x1 = interp_spec_nf[0, left_index]
            x2 = interp_spec_nf[0, left_index + 1]
            y1 = interp_spec_nf[1:, left_index]
            y2 = interp_spec_nf[1:, left_index + 1]
            x = np.log10(freq)
            interpolated_vals = (x * (y2 - y1) + x2 * y1 - x1 * y2) / (x2 - x1)
            output[ifreq] = np.sum(10.**interpolated_vals[0] * (interpolated_vals[1:, None] * map_ni), axis=0)

            output[ifreq] = hp.pixelfunc.reorder(output[ifreq], n2r=True)

            if self.unit == 'TCMB':
                conversion = 1. / K_CMB2MJysr(1., 1e9 * freq)
            elif self.unit == 'TRJ':
                conversion = 1. / K_RJ2MJysr(1., 1e9 * freq)
            else:
                conversion = 1.
            output[ifreq] *= conversion

#            output.append(result)

        if len(output) == 1:
            output = output[0]
        #else:
        #    map_data = np.row_stack(output)

        self.generated_map_freqs = freqs
        self.generated_map_data = output

        return output


    def view(self, idx=0, logged=False):
        """ View generated map using healpy's mollweide projection.

        Parameters
        ----------
        idx: int
            index of map to view. Only required if you generated maps at
            multiple frequencies.
        logged: bool
            Take the log of the data before plotting. Defaults to False.

        """

        if self.generated_map_data is None:
            raise RuntimeError("No GSM map has been generated yet. Run generate() first.")

        if self.generated_map_data.ndim == 2:
            gmap = self.generated_map_data[idx]
            freq = self.generated_map_freqs[idx]
        else:
            gmap = self.generated_map_data
            freq = self.generated_map_freqs

        if logged:
            gmap = np.log2(gmap)

        hp.mollview(gmap, coord='G',
                    title='Global Sky Model %s %s' % (str(freq), self.unit))
        plt.show()

    def write_fits(self, filename):
        """ Write out map data as FITS file.

        Parameters
        ----------
        filename: str
            file name for output FITS file
        """
        hp.write_map(filename, self.generated_map_data, column_units=self.unit)

class GSMObserver2016(ephem.Observer):
    """ Observer of the Global Sky Model.

    Generates the Observed sky, for a given point on Earth.
    Applies the necessary rotations and coordinate transformations
    so that the observed 'sky' can be returned, instead of the
    full galaxy-centered GSM.

    This class is based on pyephem's Observer(). The GSM bit can be thought of
    as an 'add on' to ephem.Observer, adding the methods generate() and  view(),
    which allows all-sky images for a given point on earth to be produced.
    """

    def __init__(self):
        """ Initialize the Observer object.

        Calls ephem.Observer.__init__ function and adds on gsm
        """
        super(GSMObserver2016, self).__init__()
        self.gsm = GlobalSkyModel2016(freq_unit='MHz')
        self.observed_sky = None

        # Generate mapping from pix <-> angles
        self.gsm.generate(1000)
        self._n_pix  = hp.get_map_size(self.gsm.generated_map_data)
        self._n_side = hp.npix2nside(self._n_pix)
        self._theta, self._phi = hp.pix2ang(self._n_side, np.arange(self._n_pix))

    def generate(self, freq):
        """ Generate the observed sky for the observer, based on the GSM.

        Parameters
        ----------
        freq: float
            Frequency of map to generate, in units of MHz (default).

        Returns
        -------
        observed_sky: np.array
            Numpy array representing the healpix image, centered on zenith,
            with below the horizon masked.
        """
        self.gsm.generate(freq)
        sky = self.gsm.generated_map_data

        # Get RA and DEC of zenith
        ra_rad, dec_rad = self.radec_of(0, np.pi/2)
        ra_deg  = ra_rad / np.pi * 180
        dec_deg = dec_rad / np.pi * 180

        # Apply rotation
        hrot = hp.Rotator(rot=[ra_deg, dec_deg], coord=['G', 'C'], inv=True)
        g0, g1 = hrot(self._theta, self._phi)
        n_side = hp.npix2nside(self.gsm.generated_map_data.shape[0])
        pix0 = hp.ang2pix(n_side, g0, g1)
        sky_rotated = sky[pix0]

        # Generate a mask for below horizon
        mask1 = self._phi + np.pi / 2 > 2 * np.pi
        mask2 = self._phi < np.pi / 2
        mask = np.invert(np.logical_or(mask1, mask2))

        self.observed_sky = hp.ma(sky_rotated)
        self.observed_sky.mask = mask

        return self.observed_sky


    def view(self, logged=False, show=False, **kwargs):
        """ View the local sky, in orthographic projection.

        Parameters
        ----------
        logged: bool
            Default False, return the log2 image
        """
        sky = self.observed_sky
        if logged:
            sky = np.log2(sky)

        hp.orthview(sky, half_sky=True, **kwargs)

        if show:
            plt.show()

        return sky

    def view_observed_gsm(self, logged=False, show=False, **kwargs):
        """ View the GSM (Mollweide), with below-horizon area masked. """
        sky = self.observed_sky
        if logged:
            sky = np.log2(sky)

        # Get RA and DEC of zenith
        ra_rad, dec_rad = self.radec_of(0, np.pi / 2)
        ra_deg  = ra_rad / np.pi * 180
        dec_deg = dec_rad / np.pi * 180

        # Apply rotation
        derotate = hp.Rotator(rot=[ra_deg, dec_deg])
        g0, g1 = derotate(self._theta, self._phi)
        pix0 = hp.ang2pix(self._n_side, g0, g1)
        sky = sky[pix0]

        coordrotate = hp.Rotator(coord=['C', 'G'], inv=True)
        g0, g1 = coordrotate(self._theta, self._phi)
        pix0 = hp.ang2pix(self._n_side, g0, g1)
        sky = sky[pix0]

        hp.mollview(sky, coord='G', **kwargs)

        if show:
            plt.show()

        return sky


