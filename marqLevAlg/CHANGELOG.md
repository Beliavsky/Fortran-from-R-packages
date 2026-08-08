# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation.
- Ported the Marquardt-Levenberg optimization loop, numerical derivative
  routines, RDM convergence, partial Hessian handling, line search, and
  mixed-model example likelihood/gradient.
- Replaced legacy packed Cholesky/inverse implementations with explicit modern
  SPD Cholesky procedures while preserving packed upper-triangle interfaces.
- Added strict regression tests and examples.
