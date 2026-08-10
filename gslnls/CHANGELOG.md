# Changelog

## 0.1.0

Initial modern Fortran/FPM translation of gslnls 1.4.2.

- Added six trust-region nonlinear least-squares methods.
- Added robust IRLS losses and weighting.
- Added parameter bounds and general weight matrices.
- Added quasi-random multi-start fitting.
- Added dense and matrix-free large-system Steihaug-CG paths.
- Added numerical Jacobian/directional-second-derivative helpers.
- Added covariance, leverage, Cook-distance, log-likelihood, and confidence
  interval utilities.
- Removed runtime R/GSL/Matrix dependencies.
- Added strict GNU Fortran regression scripts, six tests, and two examples.
- User callback invocation sites use explicit module-level procedure interfaces
  for portability with `-Werror=implicit-interface`.
