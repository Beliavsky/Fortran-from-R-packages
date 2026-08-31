# Validation

The retained test suite exercises both modern wrappers and original translated native numerical kernels.

## `test_core`

- CSR triplet conversion with duplicate summation and ignored out-of-range triplets
- sorted/valid CSR invariants
- circulant wraparound indexing
- exact RW1 and RW2 precision matrices
- spherical covariance reference values
- full nearest-distance structure including stored zero diagonal distances
- exponential covariance transform on a sparse distance matrix

## `test_cholesky_eigen`

- exact Ng-Peyton sparse Cholesky on a nontrivial SPD matrix
- sparse solve residual
- log determinant against LAPACK LU determinant
- bundled symmetric ARPACK convergence and eigenpair residual

## `test_random_mle`

- Monte Carlo mean/covariance for covariance-form Gaussian simulation
- precision-form simulation checked against `Q Cov(X) = I`
- exact linear-constraint satisfaction for precision-form sampling
- Gaussian MLE with intercept and nugget-only covariance; recovers the analytic MLE mean and variance

## `test_utils`

- Euclidean, maximum, Minkowski and great-circle distance calculations
- permutation semantics
- block `gmult`
- seasonal precision construction

The clean validation build uses GNU Fortran 14.2, `-std=f2018`, runtime checks for modern modules/tests, and `-Werror=implicit-interface` for modern public sources. The mechanically converted legacy numerical files emit obsolescent labeled-DO warnings under Fortran 2018; they compile successfully without source-algorithm changes.
