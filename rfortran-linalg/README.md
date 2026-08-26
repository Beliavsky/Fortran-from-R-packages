# rfortran-linalg

`rfortran-linalg` is the optional shared linear-algebra layer for the Fortran
translations in this repository. It provides checked double-precision APIs for
square linear systems, symmetric eigendecomposition, Cholesky factors, and
symmetric-positive-definite inverses and log determinants.

The implementation uses a pinned revision of the BSD-3-Clause
[`fortran-lapack`](https://github.com/perazz/fortran-lapack) package. It does
not require system `-llapack` or `-lblas` libraries.

The lightweight `rfortran-core` package deliberately remains independent of
BLAS and LAPACK. Basic vector operations already provided by Fortran, including
the intrinsic `norm2`, are not duplicated here.
