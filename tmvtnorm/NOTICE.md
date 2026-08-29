# Attribution and provenance

This project is a modern Fortran translation of the computational core of the R package **tmvtnorm 1.7**.

Upstream package authors:

- Stefan Wilhelm
- Manjunath B G

Upstream citation:

Stefan Wilhelm and Manjunath B G (2025), *tmvtnorm: Truncated Multivariate Normal and Student t Distribution*, R package version 1.7.

Important algorithmic references identified in the upstream source include:

- Kotecha, J. H. and Djuric, P. M. (1999), Gibbs sampling for truncated multivariate Gaussian random variables.
- Geweke (1991), simulation from multivariate normal and Student-t distributions subject to linear constraints.
- Tallis (1961), moment generating function of the truncated multinormal.
- Lee (1979, 1983), moments and estimators for truncated multivariate normal/Tobit models.
- Leppard and Tallis (1989), evaluation of truncated multinormal mean/covariance.
- Manjunath B G and Stefan Wilhelm, moments of the doubly truncated multivariate normal distribution.
- Cartinhour (1990), one-dimensional marginal densities of a truncated multivariate normal.
- Johnson and Kotz (1972), block formulas for partially truncated multivariate normals.

The attached `mvtnorm-fortran` project is redistributed as a separate vendored FPM dependency and retains its own authorship, notices and GPL-2.0-only license.

The exact upstream computational R/Fortran files used for this translation are copied under `upstream/`.
