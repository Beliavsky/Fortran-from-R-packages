# Provenance

The computational definitions in this translation were derived from `roll` 1.2.1 by Jason Foster. Reference copies of the upstream R wrapper, C++ implementation, key headers, metadata, and tests are retained under `upstream/` for auditability. They are not compiled by FPM.

The upstream implementation uses Rcpp/RcppArmadillo/RcppParallel and separate online/offline worker classes. The Fortran translation expresses the same window definitions directly with deterministic free-form Fortran loops. The `online` argument is accepted for interface parity, but both values select the same batch computation. This preserves target statistics but does not reproduce parallel scheduling or floating-point operation order of every upstream worker.

`rfortran-core` is reused for the common real kind `dp` via `r_kinds`. `rfortran-linalg` is reused for linear-system solves in `roll_lm`. No copied BLAS/LAPACK source is included.

R-specific class preservation, xts/zoo attributes, names/dimnames, `.Call` wrappers, warnings, and R list construction are intentionally omitted. Fortran arrays carry the numerical results directly. Missing real observations are represented by IEEE quiet NaNs. Logical missing values use the public integer sentinel `roll_na_logical`.
