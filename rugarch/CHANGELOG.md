# Changelog

## 0.4.0-experimental - 2026-07-22

- Added mean external regressors, two-step ARCH-in-mean, variance regressors,
  and variance targeting through the extended GARCH workflow.
- Added automatic numerical Hessian, classical covariance, per-observation
  scores, Newey-West sandwich covariance, and standard errors for extended fits.
- Added conditional raw, kernel, and semi-parametric bootstrap forecast paths,
  with partial and full simulation/refit modes.
- Added serial ARFIMA rolling, bootstrap forecast, parameter-distribution,
  multi-series, and cross-validation workflows.
- Added distribution skewness/kurtosis, GH transformation, unconditional mean,
  volatility half-life, forecast performance, and numeric index utilities.
- Added a seventh strict runtime-check test for the completed numerical surface.
- Parallel rolling, resumable rolling, R S4/S3 infrastructure, and plotting
  remain intentionally omitted.

## 0.3.0 - 2026-07-23

- Added numerical Hessian, classical covariance, Newey-West covariance, and
  sandwich covariance utilities.
- Added information criteria, weighted portmanteau, ARCH-LM, Nyblom,
  sign-bias, and adjusted Pearson goodness-of-fit diagnostics.
- Added VaR duration, GMM, Hong-Li, and Model Confidence Set tests.
- Added stationary and circular fixed-block bootstrap primitives.
- Added sequential multi-series fitting/forecasting, rolling forecasts, and
  parametric simulation/refit distributions.
- Added ARFIMA forecasts and approximate automatic AIC/BIC order selection.
- Added the complete family of rugarch forecast loss functions.
- Added a strict runtime-check test covering the new numerical surface.

## 0.2.0-experimental - 2026-07-22

- Added standardized generalized hyperbolic, normal-inverse-Gaussian, and
  generalized-hyperbolic skew-t density, CDF, quantile, and random generation.
- Added GH/NIG/GHST support to generic distribution dispatch, distribution
  fitting, GARCH likelihoods, simulation, VaR, and expected shortfall.
- Replaced preliminary FIGARCH(1,1) logic with generalized p/q polynomial
  weights and direct likelihood fitting.
- Added direct component-GARCH and arbitrary-order realGARCH fitting.
- Added joint realGARCH return/measurement likelihood and simulation.
- Added Hentschel fGARCH fitting for GARCH, TGARCH, AVGARCH, NGARCH, NAGARCH,
  APARCH, ALLGARCH, and GJR-GARCH restrictions.
- Added strict runtime-check tests for the new distributions and model fitters.

## 0.1.0-experimental - 2026-07-22

- Initial standalone modern Fortran/fpm translation.
- Added major GARCH recursions, selected fitters, distributions, ARFIMA
  utilities, forecasting, risk measures, and forecast backtests.
- Added smoke tests, examples, coverage map, and provenance documentation.
