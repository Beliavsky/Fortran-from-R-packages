# Changelog

## 0.2.0 - 2026-07-23

- Generalized DCC fitting to arbitrary scalar `(p,q,g)` orders.
- Added multivariate Normal, standardized Student-t, and Laplace densities and
  simulation, including Bessel-K support for the Laplace likelihood.
- Added Student shape estimation and Laplace DCC estimation/simulation.
- Added grouped FDCC fitting, forecasting, and simulation.
- Added static/dynamic Gaussian and Student copula fitting, transformations,
  and simulation.
- Added a Gaussian-margin Copula-GARCH fit/filter/simulation workflow.
- Added RADICAL ICA and a fitted GO-GARCH(1,1) workflow with filtering,
  forecasting, simulation, higher moments, and rolling forecasts.
- Added weighted-margin paths, FFT grid convolution/density/CDF/quantiles,
  conditional scenarios, and rolling raw-return DCC covariance forecasts.
- Added robust Huber VARX fitting and expanded diagnostic coverage.
- Added an extended runtime-checked regression suite and automated source
  license-header verification.

## 0.1.0 - 2026-07-23

- Initial experimental modern Fortran release.
- Added marginal Gaussian GARCH(1,1), DCC/ADCC, initial FDCC/copula routines,
  FastICA, GO-GARCH moment utilities, VARX, and diagnostics.
- Added GPL-3.0-only license file and headers to all Fortran sources.
