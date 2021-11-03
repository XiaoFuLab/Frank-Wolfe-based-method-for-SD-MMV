FGNSR -- A Fast Gradient method for Nonnegative Sparse Regression
=================================================================

Given a nonnegative m-by-n matrix M, FGNSR computes an n-by-n matrix X
with entries in [0,1] such that

  M * X \approx M,

and with the property that X may have only a few rows of large norm.  For
example, FGNSR can be used to compute near-separable, nonnegative matrix
factorizations.  A full description with application in hyperspectral imaging
(endmember detection) can be found here:

Gillis, N., & Luce, R. (2018). A Fast Gradient Method for Nonnegative Sparse
Regression With Self-Dictionary. IEEE Transactions on Image Processing, 27(1),
24–37.

http://doi.org/10.1109/TIP.2017.2753400


Installation
------------

FGNSR is written in Matlab and C.  A mex compatible compiler must be installed
on your machine.  To get started, run from within Matlab the script
'setup_fgnsr', located in the root directory.
