# fastmatrix-fortran

Modern Fortran translation of the numerical core of the R package `fastmatrix` 0.6-6.

## Implemented numerical areas

- Structured matrices: AR(1), compound-symmetry, circulant, Hankel, Frank, Helmert, commutation, duplication and symmetrizer matrices.
- Matrix products and transforms: `vec`, `vech`, Hadamard, Kronecker, bracket and symmetric products, Householder matrices.
- Factorizations and linear algebra: LU, LU solve/inverse, LDL, Jacobi symmetric eigensolver, matrix square root, whitening, modified Cholesky, Cholesky rank-one update/downdate, Sherman-Morrison, rank-one update, sweep operator and equilibration.
- Iterative algorithms: conjugate gradient, Gauss-Seidel, power method and Krylov matrices.
- Matrix utilities: Frobenius/1/infinity norms, condition estimate, matrix polynomial, symmetric matrix exponential/logarithm/power, batched rank-3 matrix multiplication.
- Statistics: central moments, geometric mean, skewness, kurtosis, weighted covariance, MSSD covariance, median centering, concordance correlation coefficient, Mahalanobis distances, chi d/p/q/r and Jarque-Bera statistic.
- Regression: OLS, conjugate-gradient OLS, LAPACK QR OLS, LAPACK SVD OLS and ridge regression.
- Random generation: multivariate normal, uniform sphere and uniform ball.
- Statistical compatibility: Mardia multivariate skewness/kurtosis and normality test; Harris variance-homogeneity tests (Wald, log, robust, log-robust).
- LAPACK compatibility: real nonsymmetric Schur decomposition and thin SVD.
- Matrix functions: callback-based Parlett recurrence for real upper-triangular matrices, matching upstream `matrix.fun`.
- Miscellaneous algorithms: Floyd-Warshall shortest paths and Bezier curves via de Casteljau.

The public umbrella module is `fastmatrix`.

```fortran
use fastmatrix
```

## Build

With Fortran Package Manager:

```text
fpm test
fpm run --example demo
```

The source is free-form Fortran 2018 and uses the public `dp` kind. LAPACK and BLAS are linked for Schur, SVD, and QR compatibility routines.

## Scope

R formula/data-frame/S3 presentation infrastructure is intentionally excluded. The remaining omissions are predominantly R formula/model-frame/S3 presentation infrastructure and small metadata/convenience wrappers. The numerical Schur, QR/SVD OLS, Mardia, Harris, and triangular `matrix.fun` targets are implemented in v0.2.0.
