PyGSM
=====

`PyGSM` is a Python interface for the Global Sky Model (GSM) of [Oliveira-Costa et. al., (2008)](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2966.2008.13376.x/abstract).
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

You may then install this by running `python setup.py install`.

To get a quick feel of what `PyGSM` does, have a look at the 
[quickstart guide](http://nbviewer.ipython.org/github/telegraphic/PyGSM/blob/master/docs/pygsm_quickstart.ipynb).

Q & A
-----

**Q. What's the difference between this and the `gsm.f` from the main GSM website?**
     The `gsm.f` is a very basic Fortran code, which reads and writes values to and from
     ASCII files, and uses a command line interface for input. If you want to run this code
     on an ancient computer with nothing by Fortran installed, then `gsm.f` is the way to go. 
     In contrast, `PyGSM` is a Python code that leverages a lot of other Packages so that you 
     can do more stuff more efficiently. For example: you can view a sky model in a healpy 
     image; you can write a sky model to a Healpix FITS file; and believe it or not, the 
     Python implementation is *much faster*. Have a look at the 
     [quickstart guide](http://nbviewer.ipython.org/github/telegraphic/PyGSM/blob/master/docs/pygsm_quickstart.ipynb)
     to get a feel for what `PyGSM` does.

**Q. Are the outputs of `gsm.f` and `pygsm` identical?** At the moment: **no**. The cubic
     spline interpolation implementation differs, so values will differ by as much as 
     a few percent. The interpolation code used in `gsm.f` does not have an open-source
     license (it's from [Numerical Recipes](http://www.nr.com/licenses/) ), so we haven't 
     implemented it (one could probably come up with an equivalent that didn't infringe).
     Nevertheless, the underlying PCA data are identical, and I've run tests to check that
     the two outputs are indeed comparable. 

**Q. Why is this package so large (~150 MB)?**
     The package size is dominated by the PCA healpix maps, which have about 3 million points each.
     They're compressed using HDF5 LZF, so are actually about 3x smaller than the `*.dat`
     files that come in the original `gsm.tar.gz` file. The next biggest thing is test data,
     so that the output can be compared against precomputed output from `gsm.f`.

**Q. Why do I need h5py?**
     `h5py` is required to read the PCA data, which are stored in a HDF5 file. Reading from
     HDF5 into Python is incredibly efficient, and the compression is transparent to the end user.
     This means that you can't eyeball the data using `vim` or `less` or a text editor, but if
     you're trying to do that on a file with millions of data points you're doing science wrong anyway.
   

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
