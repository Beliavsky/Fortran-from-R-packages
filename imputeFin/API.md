# API

All real data use `real(dp)`, where `dp = kind(1.0d0)`. Missing values are IEEE
quiet NaNs.

## Options and results

- `type(ar1_options)`: AR restrictions, outlier handling, tolerances, and
  Student-t SAEM controls.
- `type(imputation_options)`: sample count, burn-in, thinning, rolling window,
  and deterministic RNG seed.
- `type(ar1_fit_result)`: `phi0`, `phi1`, `sigma2`, `nu`, convergence status,
  missing/outlier indices, optional iterates, and optional conditional moments.
- `type(imputation_result)`: `values(time, series, sample)` and fitted models.
- `type(var_t_options)`: VAR order, missing-row policy, tolerance, iterations.
- `type(var_t_result)`: `phi0`, `phi(:,:,lag)`, scatter, `nu`, completed data,
  convergence information, and number of regression rows used.

## Procedures

### `fit_ar1_gaussian(y, result [, options, return_iterates, return_conditional])`

Fits a Gaussian AR(1) model to a vector or each column of a matrix. Supports
fixed random-walk coefficient (`random_walk`) and fixed zero intercept
(`zero_mean`). Missing-data estimation uses EM.

### `impute_ar1_gaussian(y, result [, options, impute_options_in])`

Draws one or more Gaussian conditional imputations. Matrix input is imputed
column by column, matching the univariate behavior of the R package.

### `impute_rolling_ar1_gaussian(y, y_imputed [, options, rolling_window, seed, status])`

Imputes a vector or matrix in rolling blocks.

### `fit_ar1_t` and `impute_ar1_t`

Student-t AR(1) estimation and latent-scale Gibbs imputation. Set
`fast_and_heuristic=.false.` to use the Gibbs/SAEM estimation path.

### `fit_var_t(y, result [, options])`

Fits a Student-t VAR(p). `phi(:,:,1)` is the first lag matrix. If
`omit_missing=.true.`, rows whose response or lags contain missing values are
excluded. Otherwise missing values are initialized by column means and updated
using conditional contemporaneous means.

### `impute_ohlc`

Imputes log close, open-minus-close, high-minus-close, and low-minus-close in
that order. The returned high is constrained not below close and the returned
low is constrained not above close.

### `impute_vol`

Log-transforms, imputes, and exponentiates a positive volume series.

### `is_inner_na` and `any_inner_na`

Return the mask or presence of NaNs located between the first and last observed
values.
