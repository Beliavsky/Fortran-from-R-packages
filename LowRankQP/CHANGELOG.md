# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of LowRankQP 1.0.6.
- Translated primal-dual predictor/corrector interior-point method.
- Translated LU, Cholesky, SMW, and PFCF solution paths.
- Replaced BLAS/LAPACK dependency with standalone Fortran linear algebra.
- Added typed options/results, objective helper, tests, examples, and strict
  GNU Fortran validation scripts.
