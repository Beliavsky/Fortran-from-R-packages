# tsDyn Modern Fortran

A modern Fortran numerical translation of substantial computational portions of the R package `tsDyn` 11.0.5.2.

The project deliberately excludes plotting, S3 method dispatch, R formulas, `ts`/`zoo` metadata, data-frame conversion, printing, tidying, LaTeX generation, GUIs, and parallel R orchestration. Numerical coverage and remaining omissions are listed below and in `API_MAP.md`.

## Implemented model families

- Univariate AR models in levels, differences, and ADF form
- Two- and three-regime SETAR models with TAR or MTAR transitions
- Two-regime additive LSTAR models
- Local-linear autoregression with direct, box-index, or automatic radius-neighbor search
- Multivariate VAR models in levels, differences, and ADF form
- Rank-restricted VECM models using fixed cointegration vectors, Engle-Granger two-step estimation, or Johansen-style maximum likelihood
- Two- and three-regime TVAR models with arbitrary transition weights
- Rank-one TVECM models with fixed or two-variable grid-searched cointegration vectors

## Implemented numerical workflows

- Simulation, filtering/fitted values, residuals, likelihood criteria, and deterministic forecasts
- AR/VAR lag selection, SETAR order selection, and VECM rank selection
- Characteristic roots and VAR/VECM coefficient representations
- Orthogonal and non-orthogonal impulse responses and VAR forecast-error variance decomposition
- Regime-specific IRFs and simulation-based GIRFs for SETAR, TVAR, and TVECM
- Residual, wild, IID, circular-block, Rademacher, Mammen, and Gaussian-multiplier resampling
- AR and VAR bootstrap forecast paths
- Expanding- or fixed-window rolling forecasts for AR, SETAR, and VAR
- Correlation-integral delta statistics and shuffle/linear-surrogate tests
- BBC (2004)-style threshold unit-root LR/Wald/LM statistics
- Kapetanios-Shin (2006)-style sup, average, and exponential-average statistics
- Johansen trace and maximum-eigenvalue statistics from fitted ML VECMs
- Forecast error measures: ME, RMSE, MAE, MPE, and MAPE

## Build and test

Requirements:

- GNU Fortran or another Fortran 2018 compiler
- LAPACK
- BLAS

Runtime-checked build and tests:

```sh
make debug-check
```

Optimized build and tests:

```sh
make release-check
```

Both targets compile with warnings treated as errors. The debug target additionally enables `-fcheck=all`.

An `fpm.toml` manifest is included. `fpm` was not installed in the validation environment, so the manifest is not claimed as tested.

## Applications

```sh
build/debug/bin/demo_tsdyn
build/debug/bin/fit_csv data/sample.csv ar 2
build/debug/bin/fit_csv data/sample.csv var 1
build/debug/bin/threshold_example
```

`fit_csv` expects a header and a first date/text column followed by one or more numeric series. It supports `ar`, `setar`, `lstar`, `var`, and fixed-vector `vecm` modes.

## Important numerical differences

- LAPACK least squares and eigensolvers replace R's linear-algebra wrappers.
- LSTAR fitting uses deterministic grid search followed by local pattern refinement, not R's optimizer stack.
- Johansen estimation is implemented directly from residual covariance matrices; normalization and finite-precision endpoints can differ from `urca`.
- Local-linear autoregression uses direct or dynamic box-index radius searches rather than the original C storage layout.
- Bootstrap and simulation use the Fortran intrinsic random-number generator, not R's RNG streams.
- Threshold tests return statistics and selected thresholds but do not reproduce every package bootstrap critical-value procedure.

Exact R output, optimizer endpoints, random streams, and external-package behavior are not claimed.

## Explicit omissions

The following computational families remain unimplemented:

- Additive nonlinear AR fitting through `mgcv::gam` (`aar`)
- Neural-network time-series fitting and selection through `nnet` (`nnetTs`, `selectNNET`)
- The generic multi-transition STAR construction and gradient/parameter-management engine beyond the tested two-regime LSTAR implementation
- Full Hansen (1999), Hansen-Seo (2002), and Seo (2006) bootstrap critical-value algorithms
- Symbolic VECM expressions
- Full analytical coefficient covariance/confidence-interval methods for every model
- Nonlinear FEVD beyond VAR FEVD and regime/GIRF response calculations

These omissions are not represented as implemented features.

## License

The original package declares `GPL (>= 2)`. This translation is licensed as `GPL-2.0-or-later`. The full GNU GPL version 2 text is in `LICENSE`, and every Fortran source file contains a machine-checked SPDX and license notice.
