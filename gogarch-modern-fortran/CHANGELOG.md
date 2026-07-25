# Changelog

## 0.2.0 - 2026-07-23

- Added GARCH(p,q) filtering, fitting, forecasting, and simulation.
- Added APARCH(p,o,q) powered-scale recursions and leverage estimation.
- Added standardized Normal, skew-Normal, Student-t, skew-Student, GED, and
  skew-GED densities and random generation.
- Added conditional shape, skew, and APARCH delta estimation.
- Propagated factor specifications through ICA, MM, NLS, ML, direct-angle
  construction, likelihood evaluation, forecasting, and fitted simulation.
- Added full higher-order factor parameter extraction.
- Extended the CSV application with model, distribution, and order arguments.
- Added two numerical test suites for the new univariate and GO-GARCH paths.

## 0.1.0 - 2026-07-23

- Initial modern Fortran computational translation.
- Added FastICA, methods-of-moments, nonlinear least-squares, and
  Euler-angle maximum-likelihood GO-GARCH estimation.
- Added Gaussian GARCH(1,1) factor models, covariance paths, forecasts,
  fitted simulation, tests, examples, and license checks.
