# API Reference

All real-valued numerical APIs use `real(dp)`, where
`dp = kind(1.0d0)`. Procedures return an explicit status code when failure is
possible.

## Main result and option types

### `type(acd_order)`

- `p`: duration/residual lag order
- `r`: innovation/residual moving-average order used by AMACD/TAMACD
- `q`: conditional-mean lag order

### `type(acd_fit_options)`

Controls model code, distribution code, lag order, forced unit-mean errors,
iteration limits, restarts, tolerance, RNG seed, Hessian calculation, and robust
standard errors.

### `type(acd_fit_result)`

Contains:

- all parameters
- model and distribution parameter subvectors
- conditional means and residuals
- log likelihood, AIC, BIC, and MSE
- Hessian and conventional covariance matrix
- conventional standard errors
- robust sandwich covariance and standard errors when requested
- convergence, iteration, evaluation, and status information

## Model procedures

### `filter_acd`

Filters a positive duration vector under one of the 13 model codes. Optional
arguments support thresholds, spline breakpoints, exogenous regressors, and daily
recursion resets.

### `acd_loglik`

Returns the complete conditional log likelihood and the corresponding conditional
means and multiplicative residuals.

### `simulate_acd`

Simulates any model/distribution combination. It accepts optional burn-in,
user-supplied innovations, starting durations and means, breakpoints, and
exogenous regressors.

### Model constants

`MODEL_ACD`, `MODEL_LACD1`, `MODEL_LACD2`, `MODEL_EXACD`, `MODEL_AMACD`,
`MODEL_ABACD`, `MODEL_AACD`, `MODEL_TACD`, `MODEL_BACD`, `MODEL_BCACD`,
`MODEL_SNIACD`, `MODEL_LSNIACD`, and `MODEL_TAMACD`.

`model_code`, `model_name`, `model_parameter_count`, and
`default_model_parameters` provide metadata and defaults.

## Estimation and inference

### `acd_fit_model`

Fits a model by bounded Nelder-Mead maximum likelihood with optional random
restarts. Optional vectors provide starting values, lower/upper bounds, and a
logical fixed-parameter mask.

### `acd_score_matrix`

Computes observation-level numerical likelihood scores.

### `acd_bootstrap_se`

Performs a residual bootstrap, refits each simulated series, and returns
bootstrap standard errors.

### `forecast_acd`

Produces deterministic conditional-mean forecasts by propagating unit-mean
future innovations through the fitted recursion.

### `default_parameter_bounds`

Returns broad, model-aware numerical bounds suitable for the built-in optimizer.
Users can replace any bounds when fitting.

## Distribution API

### Generic procedures

- `distribution_parameter_count`
- `distribution_pdf`
- `distribution_logpdf`
- `distribution_cdf`
- `distribution_quantile`
- `sample_distribution`
- `forced_scale`

The generic functions accept a distribution code, a parameter vector, and a
logical `force_mean` flag.

### Distribution constants

`DIST_EXPONENTIAL`, `DIST_WEIBULL`, `DIST_BURR`, `DIST_GENGAMMA`,
`DIST_GENF`, `DIST_QWEIBULL`, `DIST_MIXQWE`, `DIST_MIXQWW`,
`DIST_MIXINVGAUSS`, and `DIST_BIRNBAUM_SAUNDERS`.

### Named distribution procedures

Burr:
`dburr`, `pburr`, `qburr`, `rburr`, `burr_expectation`.

Generalized gamma:
`dgengamma`, `pgengamma`, `qgengamma`, `rgengamma`, `gengamma_hazard`.

Generalized F:
`dgenf`, `pgenf`, `qgenf`, `rgenf`, `genf_hazard`.

q-Weibull:
`dqweibull`, `pqweibull`, `qqweibull`, `rqweibull`,
`qweibull_expectation`, `qweibull_hazard`.

Mixtures:
`dmixqwe`, `pmixqwe`, `qmixqwe`, `rmixqwe`, `mixqwe_hazard`,
`dmixqww`, `pmixqww`, `qmixqww`, `rmixqww`, `mixqww_hazard`,
`dmixinvgauss`, `pmixinvgauss`, `qmixinvgauss`, `rmixinvgauss`, and
`mixinvgauss_hazard`.

Birnbaum-Saunders:
`dbirnbaum_saunders`, `pbirnbaum_saunders`,
`qbirnbaum_saunders`, and `rbirnbaum_saunders`.

## Transaction and diurnal APIs

### `compute_durations`

Constructs trade, price, or volume durations from parallel arrays containing
calendar dates and seconds since midnight. The result includes event timestamps,
durations, prices, aggregate volume, and transaction counts.

Duration constants:
`DURATION_TRADE`, `DURATION_PRICE`, and `DURATION_VOLUME`.

### `diurnal_adjust`

Fits a time-of-day scale and returns fitted scales, adjusted durations, and a
regular diagnostic grid. Group identifiers can request separate fits by weekday,
date, or any user-defined grouping.

Method constants:
`DIURNAL_CUBIC_SPLINE`, `DIURNAL_SMOOTH_SPLINE`,
`DIURNAL_SUPER_SMOOTHER`, and `DIURNAL_FFF`.

Low-level smoothing procedures `flexible_fourier_fit` and
`super_smoother_fit` are also public.

## Diagnostics

- `standardize_residuals`: probability-integral or Cox-Snell transform
- `acf_acd`: autocorrelation coordinates and confidence limit
- `residual_density_acd`: Gaussian KDE and fitted theoretical density
- `qqplot_acd`: empirical and theoretical quantile coordinates
- `summarize_durations`: summary statistics and rolling mean
- `rolling_mean`: efficient moving mean
- `hazard_diagnostics`: nonparametric and model-implied hazard coordinates
- `likelihood_profile`: one- or two-parameter log-likelihood grid

## Specification tests

- `test_rm_acd`: remaining-ACD LM test
- `test_st_acd`: smooth-transition ACD alternative
- `test_tv_acd`: time-varying ACD alternative

Each test supports robust and nonrobust variants and returns
`type(lm_test_result)` with statistic, degrees of freedom, p-value, and status.

## Status codes

- `ACDM_SUCCESS`
- `ACDM_BAD_INPUT`
- `ACDM_BAD_PARAMETER`
- `ACDM_NUMERIC_FAILURE`
- `ACDM_NOT_CONVERGED`
- `ACDM_SINGULAR`
