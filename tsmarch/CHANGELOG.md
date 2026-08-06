# Changelog

## 0.1.0

- Initial modern Fortran/FPM computational translation of `tsmarch` 1.0.3.
- Added DCC, asymmetric DCC, constant-correlation, and copula-GARCH engines.
- Added FastICA, RADICAL-style ICA, GO-GARCH estimation, higher moments, and
  forecasts.
- Added covariance utilities, multivariate distributions, ESCC diagnostics,
  risk measures, and convolution distribution functions.
- Vendored compatible `ghyp`, `tsdistributions`, and `tsgarch` sources.
- Retained the user-supplied LGPL-3 `nloptr` translation as unlinked provenance
  because it is incompatible with GPL-2-only static combination.
- Added checked and optimized GNU Fortran build scripts, five tests, and a
  demonstration program.
