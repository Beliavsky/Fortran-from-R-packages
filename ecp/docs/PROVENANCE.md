# Provenance

Upstream package: **ecp** version 3.1.6, published on CRAN in August 2024.

Upstream authors named in `DESCRIPTION`:

- Nicholas A. James
- Wenyu Zhang
- David S. Matteson

Upstream license: `GPL (>= 2)`.

The `upstream/` directory retains the package metadata, R source files, selected manual pages, and original native sources used to understand and validate the translation. The maintained Fortran in `src/` is a new implementation and does not compile or vendor Rcpp/C++ into the Fortran package.

Algorithm provenance includes the methods cited by upstream for energy-based divisive/agglomerative change-point analysis, CP3O pruning/dynamic programming, and kernel change-point analysis. Function-level source paths are recorded in `fpm.toml`.

No translated dependency source, BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, or `rfortran-compat` source is copied into this package.
