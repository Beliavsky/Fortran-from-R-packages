# Changelog

## 0.3.0

Second parity-focused release.

- Added joint numeric ARIMAX/xreg fitting, forecasting, simulation support, and fixed-structure `arima_refit`.
- Added an integrated ARIMA diffuse state-space likelihood, including a missing-observation Kalman path and a regression test against the closed-form random-walk likelihood.
- Corrected seasonal/nonseasonal AR polynomial multiplication in simulation and impulse calculations.
- Tightened `auto_arima`: upstream-ordered stepwise candidate/restart traversal, configurable starts/model budget/information criterion, fixed d/D controls, and approximation `truncate` selection followed by full-sample ML refit.
- Added analytic ETS class-3 forecast variance and simulation/bootstrap prediction intervals using R type-8 quantiles.
- Added linear-process bootstrap confidence intervals for tapered ACF/PACF.
- Added BATS/TBATS forced control semantics and fixed-structure `bats_refit` / `tbats_refit`; fitted models now retain seed states and training data.
- Reworked spline GCV/fitted values around the exact constant-plus-linear unpenalized null space instead of a finite large-variance approximation.
- Added `test_parity_v03` coverage for ARIMAX/refit/truncate, diffuse and missing-value likelihoods, ETS class-3/simulation intervals, tapered bootstrap intervals, type-8 quantiles, and BATS/TBATS refitting/controls.

## 0.2.0

Parity-focused follow-up release.

- Added Gaussian stationary innovations ML fitting for ARIMA, with CSS initialization,
  corrected process-mean parameterization, and impulse-response forecast variances.
- Added OCSB seasonal differencing and expanded unit-root/differencing controls.
- Replaced the earlier seasonal-decomposition surrogate with self-contained LOESS/STL and
  iterative MSTL machinery; added `stlf`-style callback forecasting.
- Expanded ETS initial-state optimization, automatic candidate search, forbidden-model and
  Hyndman-Akram-Archibald admissibility checks, fixed smoothing-parameter handling, and
  additive-model forecast-variance recursions. Fixed `ses_fit(alpha=...)` being ignored.
- Expanded BATS/TBATS fitting/selection, ARMA-error optimization, non-integer TBATS periods,
  harmonic selection, and state-space forecast variances.
- Added generic bagged-model and `CVar` callback APIs, calendar/business-day/Easter helpers,
  future Fourier/seasonal feature helpers, tapered ACF/PACF, and GCV spline selection.
- Added `test_parity_v02` regression coverage for the new parity work.

## 0.1.0

Initial modern Fortran translation of the computational core of `forecast` 9.0.2.

Highlights include ETS, ARIMA/auto-ARIMA, ARFIMA, nnetar, BATS/TBATS, DSHW, Theta,
Croston, benchmark forecasts, Box-Cox utilities, time-series features, accuracy/DM tests,
regression/modelAR, stochastic spline forecasting, cleaning/decomposition/bootstrap,
rolling-origin CV, bagged ETS and residual diagnostics.

The supplied `fracdiff`, `urca`, and `nnet` Fortran translations are bundled as local FPM
path dependencies.
