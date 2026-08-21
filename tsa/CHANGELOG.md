# Changelog

## 0.1.0 - 2026-08-18

- Initial modern free-format Fortran/FPM translation of TSA 1.3.1 computational
  routines.
- Reused supplied Fortran translations of `leaps` and `tseries`; retained
  supplied `locfit` and `mgcv` translations for provenance.
- Added AR/ARMA identification, spectra, EACF, Box-Cox profiling, diagnostics,
  outlier detection, QAR/GARCH/TAR simulation, SETAR MAIC/CLS fitting and
  threshold inference, ARIMA/ARIMAX CSS fitting, and bootstrap support.
- Added regression tests and explicit translation-status documentation.

## 0.2.0 - 2026-08-18

- Added TSA transfer functions, dynamic interventions, fixed/init masks,
  stationary exact ML, bootstrap parity, multiseries TAR, and missing-value
  handling.

## 0.3.0 - 2026-08-18

- Added diffuse integrated/seasonal ARIMA state-space ML with missing-data
  prediction updates and standardized innovations.
- Ported R's default Gardner/AS154 stationary covariance initializer.
- Added R-style BFGS optimization and TSA regression/transfer `parscale` rules.
- Added univariate and multivariate `spec.pgram` computational parity, including
  modified Daniell kernels, coherence, and phase.
- Added the default Yule-Walker `spec.ar` path and dedicated spectral/diffuse
  regression tests.
- Retained bounded-memory and overflow-safe fallbacks for strict IEEE checking.

## 0.4.0 - 2026-08-18

- Added R/TSA SVD rotation for multiple all-free regression columns, with
  coefficient and covariance mapping back to the original xreg basis.
- Added a direct one-sided Jacobi SVD for right singular vectors rather than
  forming `X^T X`, preserving conditioning for nearly collinear regressors.
- Ported R's `optimHess` numerical derivative-of-gradient construction and
  switched ARIMAX covariance estimation to that Hessian path.
- Replaced the numerical AR/SAR covariance-transform Jacobian with TSA's
  analytic Levinson/PACF derivative recursion.
- Added a self-contained pivoted-LU inverse for Hessian covariance solves,
  following the factor/solve structure of LAPACK general solves without adding
  a mandatory external BLAS/LAPACK dependency.
- Matched ordinary `stats::arima` differenced-regression `parscale` behavior
  while preserving TSA's specialized transfer/IO scaling branch.
- Added complete-case regression initialization for ML fits containing missing
  response or xreg values.
- Added `test_xreg_parity` and `test_linalg_parity`; extended optimizer tests
  with scaled `optimHess` checks.

## 0.5.0 - 2026-08-18

- Added non-default `spec.ar` fitting methods: Burg/Burg2, OLS, and ML, with
  fixed-order and AIC order-selection paths; the Burg recursion is a direct
  translation of R's current `stats/src/burg.c` update order.
- Added R-style compact symmetric `tskernel` coefficients and constructors for
  Daniell, modified Daniell, Fejer, and Dirichlet kernels, including composite
  vector-`m` Daniell and modified-Daniell kernels.
- Added arbitrary compact-kernel smoothing to univariate and multivariate
  `spec_pgram`.
- Added per-series taper fractions for multivariate periodograms, including
  series-specific taper-energy and degrees-of-freedom corrections while
  preserving coherence scaling.
- Added `test_spectral_methods` covering the new AR methods, AIC selection,
  kernel expansion/normalization, compact-vs-expanded smoothing, and
  per-series taper behavior.
- Added an R-side v0.5 spectral parity harness.
