# rfortran-linalg

`rfortran-linalg` is the optional shared linear-algebra layer for the Fortran
translations in this repository. It provides checked double-precision APIs for
square linear systems, symmetric eigendecomposition, Cholesky factors,
direct and prefactored symmetric-positive-definite solves, and SPD inverses
and log determinants. It also provides real and complex thin singular-value
decomposition, economy-size QR, rank-revealing pivoted-QR column selection,
general matrix inversion,
singular-values-only decomposition, SVD-based numerical rank, full-rank QR
least squares, rank-revealing SVD least squares, full SVD, general-matrix
spectral radius, general real and complex eigendecomposition, and real Schur
decomposition. The real API offers an eigenvalues-only path as well as right
eigenvectors; the latter preserves LAPACK's real encoding of conjugate pairs.
Square solves support both real and complex vectors and matrices. Real and
complex matrix balancing and complex Schur decomposition are also available.
Pivoted-QR least squares returns the numerical rank and supports vector or
matrix right-hand sides with a caller-selectable rank threshold.
Checked triangular inversion supports either triangle and unit or non-unit
diagonals. General determinants and signed log-absolute-determinants use LU
factorization and report singularity through the usual `info` result.

The implementation uses a pinned revision of the BSD-3-Clause
[`fortran-lapack`](https://github.com/perazz/fortran-lapack) package. It does
not require system `-llapack` or `-lblas` libraries.

The lightweight `rfortran-core` package deliberately remains independent of
BLAS and LAPACK. Basic vector operations already provided by Fortran, including
the intrinsic `norm2`, are not duplicated here.
