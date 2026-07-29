# API reference

All public user-facing procedures are available from:

```fortran
use bekks
```

The real kind is `dp = kind(1.0d0)`.

## Model constants

```fortran
bekk_full
bekk_diagonal
bekk_scalar
```

Status constants:

```fortran
bekk_ok
bekk_invalid_input
bekk_invalid_parameters
bekk_linalg_failure
bekk_no_convergence
```

A fit that reaches the iteration limit can return `bekk_no_convergence` while still containing usable estimates, filtered covariances, scores, and diagnostics.

## Core types

### `bekk_spec_type`

- `model_type`
- `asymmetric`
- `signs(:)`
- `initial_theta(:)`

Construct with:

```fortran
spec=bekk_spec(model_type,asymmetric,signs,initial_theta)
```

All arguments are optional.

### `bekk_parameters`

Contains the lower-triangular intercept matrix `c`, full matrices `a`, `b`, `g`, and scalar coefficients for scalar BEKK.

### `bekk_fit_result`

Contains:

- specification and unpacked parameters;
- packed parameter vector;
- OPG and robust covariance matrices;
- standard errors and t values;
- per-observation score matrix;
- filtered covariance path and standardized residuals;
- likelihood path, log likelihood, AIC, and BIC;
- stationarity, convergence, status, and iteration count.

### Forecast, risk, and diagnostic results

- `bekk_filter_result`
- `bekk_forecast_result`
- `bekk_var_result`
- `bekk_backtest_result`
- `bekk_virf_result`
- `bekk_portmanteau_result`

## Parameter packing

```fortran
integer function parameter_count(n,model_type,asymmetric)
subroutine unpack_parameters(theta,n,model_type,asymmetric,parameters,status)
function pack_parameters(parameters) result(theta)
```

Packing order:

1. lower triangle of `C`, column by column;
2. full or diagonal `A`, or scalar `a`;
3. `B` when asymmetric;
4. full or diagonal `G`, or scalar `g`.

## Starting values

```fortran
subroutine initial_parameters(data,model_type,asymmetric,parameters,status)
subroutine random_initial_parameters(data,model_type,asymmetric,state,parameters,status,n_trials)
```

Original-style packed-vector wrappers:

```fortran
call grid_search_bekk(data,theta,likelihood,status)
call grid_search_asymmetric_bekk(data,signs,theta,likelihood,status)
call grid_search_dbekk(data,theta,likelihood,status)
call grid_search_asymmetric_dbekk(data,signs,theta,likelihood,status)
call grid_search_sbekk(data,theta,likelihood,status)
call grid_search_asymmetric_sbekk(data,signs,theta,likelihood,status)
```

Each has a corresponding `random_grid_search_*` form accepting an explicit `rng_state` and optional `n_trials`.

## Estimation and inference

```fortran
call bekk_fit(spec,data,result,max_iter,criterion,use_qml)
call bhhh_fit(data,spec,result,max_iter,criterion,use_qml)
```

`use_qml=.true.` selects robust sandwich standard errors; it does not change the Gaussian QML objective.

Low-level score and Hessian interfaces exist for every model:

```fortran
call score_bekk(theta,data,score,status)
call score_asymm_bekk(theta,data,signs,score,status)
call score_dbekk(...)
call score_asymm_dbekk(...)
call score_sbekk(...)
call score_asymm_sbekk(...)

call hesse_bekk(theta,data,hessian,status)
! analogous full, diagonal, scalar, and asymmetric variants
```

General inference helper:

```fortran
call qml_covariance(theta,data,model_type,asymmetric,signs, &
  covariance,robust,status)
```

`covariance` is the inverse outer-product-of-gradients estimate. `robust` is the Hessian/OPG sandwich estimate.

## Filtering and likelihoods

Likelihood functions:

```fortran
loglike_bekk
loglike_asymm_bekk
loglike_dbekk
loglike_asymm_dbekk
loglike_sbekk
loglike_asymm_sbekk
```

Filtering wrappers:

```fortran
call sigma_bekk(theta,data,result)
call sigma_bekk_asymm(theta,data,signs,result)
call sigma_dbekk(...)
call sigma_dbekk_asymm(...)
call sigma_sbekk(...)
call sigma_sbekk_asymm(...)
```

The filter result contains `h(N,N,T)` and standardized residuals `residuals(T,N)`.

## Stationarity and covariance

```fortran
valid_bekk(theta,n)
valid_asymm_bekk(theta,n,expected_indicator)
valid_dbekk(...)
valid_asymm_dbekk(...)
valid_sbekk(...)
valid_asymm_sbekk(...)

expected_indicator_value(data,signs)
indicator_function(return_vector,signs)
call unconditional_covariance(parameters,expected_indicator,h,status)
call covariance_to_volatility(h,standard_deviation,correlation)
```

## Simulation

Convenience wrappers:

```fortran
call simulate_bekk_full(...)
call simulate_bekk_asymm(...)
call simulate_dbekk(...)
call simulate_dbekk_asymm(...)
call simulate_sbekk(...)
call simulate_sbekk_asymm(...)
```

The general interface also accepts a fixed innovation matrix:

```fortran
call simulate_bekk_model(theta,nobs,n,model_type,asymmetric,state,data,h,status, &
  signs,expected_indicator,innovations)
```

`innovations`, when supplied, has shape `nobs x n` and makes simulation deterministic.

## Forecasting

```fortran
call forecast_bekk(fit,n_ahead,result,confidence_level)
```

Returns forecast covariance matrices, standard deviations, correlations, and lower/upper covariance paths.

## Volatility impulse responses

```fortran
call virf_bekk(fit,h0,shock,periods,result,confidence_level)
```

`shock` is in standardized-innovation units. Results are returned in lower-triangle `vech` order. Confidence intervals use a numerical Jacobian and the fitted parameter covariance matrix.

## Value at risk

```fortran
call var_bekk_fit(fit,p,result,portfolio_weights,distribution)
call var_bekk_forecast(fit,forecast,p,result,portfolio_weights,distribution)
```

Supported residual distributions:

- `'normal'`
- `'empirical'`
- `'t'`

Without `portfolio_weights`, marginal VaR is returned for each series. With weights, one portfolio VaR series is returned.

The convention follows the source package: `p=0.99` requests the one-percent lower-tail VaR, which is normally negative for returns.

## Backtesting

```fortran
call coverage_tests(returns,var,p,kupiec_stat,kupiec_pvalue, &
  christoffersen_stat,christoffersen_pvalue)
call backtest_forecasts(returns,var,p,result)
call rolling_backtest(data,spec,window_length,p,n_ahead,result, &
  portfolio_weights,distribution,max_iter)
```

## Diagnostics and Monte Carlo evaluation

```fortran
call portmanteau_test(fit,lags,result)
call bekk_mc_eval(theta_true,spec,sample_sizes,iterations,state,mse,status,max_fit_iter)
rmse_parameters(theta_estimated,theta_true)
```

## Matrix and numerical utilities

```fortran
elimination_mat
commutation_mat
duplication_mat
diag_selection_mat
cut_mat_symmetric
cut_mat_asymmetric
vech_lower
unvech_lower
y_lag_cr
extract_csd
general_inverse
symmetric_sqrt
```
