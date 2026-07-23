# rmgarch modern Fortran

**Status: experimental, partial translation.** This project is an independent
modern Fortran reimplementation of computational methods from R `rmgarch`
1.4-2. It is not a drop-in replacement and does not reproduce the R class,
data-indexing, solver, or `rugarch` ecosystems.

Every Fortran source file begins with:

```fortran
! SPDX-License-Identifier: GPL-3.0-only
```

and carries a GPL version 3 only notice. The complete license text is in
`LICENSE`.

## Tested computational coverage

The checked test suite exercises the following implemented paths:

- Gaussian univariate GARCH(1,1) fitting and filtering
- Two-step raw-return DCC fitting
- General DCC(p,q,g)/ADCC filtering, constrained estimation, history-aware
  forecasting, simulation, persistence, and rolling forecasts
- Gaussian, standardized Student-t, and multivariate Laplace DCC likelihoods;
  Student shape estimation is tested
- Grouped FDCC(1,1) parameter construction, filtering, estimation,
  forecasting, and simulation
- Static and dynamic Gaussian and Student copula densities, likelihoods,
  estimation, transformations, and simulation
- A model-level Copula-GARCH workflow with Gaussian GARCH(1,1) margins,
  dynamic Gaussian copula fitting, filtering, and fitted simulation
- FastICA and pairwise-rotation RADICAL ICA
- Square GO-GARCH approximation with ICA, independent Gaussian GARCH(1,1)
  factors, filtering, forecasting, fitted simulation, and rolling forecasts
- GO-GARCH covariance, correlation, volatility, co-skewness, co-kurtosis,
  portfolio moments, and covariance/co-skewness/co-kurtosis betas
- Multivariate Normal, standardized Student-t, and multivariate Laplace
  log densities and simulation
- Weighted portfolio margins and time-varying weighted-margin paths
- Radix-2 FFT grid convolution with density, CDF, quantile, and moment queries
- Conditional DCC scenario generation, first through fourth scenario moments,
  and portfolio scenario projection
- VARX fitting, filtering, forecasting, simulation, and Huber robust fitting
- Correlation distances, DCC constancy diagnostics, and Mardia diagnostics

`API_MAP.md` maps the R-facing concepts to the actual Fortran procedures.
`VALIDATION.md` lists the exact test coverage and commands.

## Important limitations

The following are **not** claimed as implemented:

- Full `rugarch` marginal compatibility: ARFIMA/ARMA mean models, external
  regressors throughout, EGARCH/APARCH/GJR/realGARCH families, skewed marginal
  distributions, and the full univariate solver surface
- Standard errors, Hessians, inference tables, all original solver choices,
  fixed-parameter machinery, and exact R optimization behavior
- Semiparametric SPD marginals and the complete Copula-GARCH marginal
  transformation pipeline
- GO-GARCH NIG/GH component-distribution estimation, arbitrary rectangular
  factor structures, and the original specialized NIG/GH FFT inversion
- Parallel rolling orchestration, serialized refits, resume/recovery logic,
  and date-indexed forecast objects
- Full `fScenario`/`fMoments` class behavior or exact parity with every
  post-estimation method
- Plotting, R S3/S4 classes, formula parsing, `xts`/`zoo`, and datasets

The implementation should be treated as a tested numerical starting point, not
as proof of numerical equivalence with R `rmgarch`.

## Build and test

GNU Make:

```text
make clean
make check
make all fit_csv
./build/demo_rmgarch
./build/fit_csv data/sample_returns.csv
```

The project also includes `fpm.toml`, but the release validation recorded in
`VALIDATION.md` used GNU Fortran and GNU Make because `fpm` was not installed in
the validation environment.

## Basic examples

General two-step DCC from an `nobs x nassets` return matrix:

```fortran
use rmgarch

type(multivariate_garch_fit) :: fit

fit = fit_two_step_dcc_general(returns, p=1, q=1, g=0)
```

GO-GARCH fit and forecast:

```fortran
use rmgarch

type(gogarch_fit_result) :: fit
real(dp) :: factor_sigma(5, size(returns,2))
real(dp) :: covariance(size(returns,2), size(returns,2), 5)
real(dp) :: correlation(size(returns,2), size(returns,2), 5)

fit = fit_gogarch11(returns)
call forecast_gogarch11(fit, 5, factor_sigma, covariance, correlation)
```

## Array conventions

- Observations and standardized scores: `nobs x nassets`
- Dynamic covariance/correlation: `nassets x nassets x nobs`
- Rolling covariance/correlation: `nassets x nassets x horizon x nfits`
- GO-GARCH component scales: `nobs x nfactors`
- GO-GARCH mixing matrix: `nassets x nfactors`
