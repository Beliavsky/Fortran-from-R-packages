# NOTICE

This directory is a modern Fortran translation of computational code from the R package **roll** version 1.2.1.

Upstream package:

- Package: `roll`
- Author and maintainer: Jason Foster
- Upstream repository: <https://github.com/jasonjfoster/roll>
- Upstream license: GPL (>= 2)

The translated/new Fortran source is distributed under GPL-2.0-or-later, consistent with the upstream package license. The full GPL version 2 text is included in `LICENSE`.

The package reuses the sibling `rfortran-core` package for the shared `dp` real kind and `rfortran-linalg` for checked linear solves used by rolling regressions. Those dependencies are referenced through FPM path dependencies and are not copied into this package.

No Rcpp, RcppArmadillo, RcppParallel, BLAS, LAPACK, or dependency source is vendored in this translation.
