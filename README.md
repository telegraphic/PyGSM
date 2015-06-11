PyGSM
=====

PyGSM is aPython interface for the Global Sky Model (GSM) or Oliveira-Costa et. al. (2008)
The GSM generates all-sky maps in Healpix format of diffuse Galactic radio emission
from 10 MHz to 94 GHz.

This is *not* a wrapper of the original Fortran code, it is a python-based equivalent
that has some additional features and advantages.
Instead of the original ASCII DAT files that contain the principal component analysis
(PCA), data are stored in HDF5, which is more efficient.

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
