PyGSM
=====

`PyGSM` is a Python interface for the Global Sky Model (GSM) or Oliveira-Costa et. al. (2008)
The GSM generates all-sky maps in Healpix format of diffuse Galactic radio emission
from 10 MHz to 94 GHz.

This is *not* a wrapper of the original Fortran code, it is a python-based equivalent
that has some additional features and advantages, such as healpy integration for imaging.
Instead of the original ASCII DAT files that contain the principal component analysis
(PCA), data are stored in HDF5, which can be quickly read into memory, and takes up less space to boot.

Quickstart
----------

The first thing to do will be to make sure you've got the dependencies: 

* [numpy](http://www.numpy.org/)
* [scipy](http://www.scipy.org/install.html)
* [healpy](http://healpy.readthedocs.org/en/latest/)
* [h5py](http://www.h5py.org/)
* [astropy](http://www.astropy.org/)

Then clone the directory

        git clone https://github.com/telegraphic/PyGSM

and run any scripts or iPython sessions from there (installation via `setup.py` is on the todo list).

To get a quick feel of what `PyGSM` does, have a look at the 
[quickstart guide](http://nbviewer.ipython.org/github/telegraphic/PyGSM/blob/master/pygsm_quickstart.ipynb).

References
----------

The PCA data contained here is from http://space.mit.edu/~angelica/gsm/index.html.

The original paper is:

```
A. de Oliveira-Costa, M. Tegmark, B.M. Gaensler, J. Jonas, T.L. Landecker and P. Reich
A model of diffuse Galactic radio emission from 10 MHz to 100 GHz
Mon. Not. R. Astron. Soc. 388, 247-260 (2008)
doi:10.111/j.1365-2966.2008.13376.x
```

Which is published in [MNRAS](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2966.2008.13376.x/abstract)
and is also available on the [arXiv](http://arxiv.org/abs/0802.1525).
