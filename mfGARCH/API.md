# API overview

All public entities are re-exported by `use mfgarch`.

## Main types

### `mfgarch_model`

Physical model parameters and structural options:

- `mu`, `alpha`, `beta`, `gamma`
- `m`, `theta`, `w1`, `w2`
- `theta_two`, `w1_two`, `w2_two`
- `k`, `k_two`
- `asymmetric`, `unrestricted_weights`, `has_second`

For restricted beta weights, set `w1 = 1`. For `k = 1`, the sole weight is one
and shape parameters are not estimated.

### `mfgarch_fit_control`

Controls optimizer method, tolerances, iteration limits, multi-stage fitting,
tracing, and covariance calculation.

### `mfgarch_fit_result`

Contains the fitted model, log likelihood, BIC, `tau`, `g`, residuals, beta
weights, variance ratio, tau forecast, and conventional/robust/OPG covariance
matrices and standard errors.

### `mfgarch_simulation`

Contains daily and intraday returns, covariates, long- and short-run variance
components, realized variances, and 5/22-period rolling averages.

## Estimation

```fortran
call fit_mfgarch(returns, period, start_model, result, status, &
  covariate=x, period_two=period2, covariate_two=x2, &
  variance_period=variance_groups, control=control)
```

Only `returns`, `period`, `start_model`, `result`, and `status` are required.
When `start_model%k > 0`, supply `covariate`. When
`start_model%has_second` is true, also supply `period_two` and
`covariate_two`.

## Components and likelihood

- `beta_weights`
- `low_frequency_log_tau`
- `build_tau`
- `calculate_g`
- `likelihood_contributions`
- `log_likelihood`
- `variance_ratio`
- `forecast_tau`
- `forecast_garch`

## Forecasting

```fortran
call predict_variance(model, horizons, tau_forecast, last_return, &
  conditional_g, conditional_tau, forecasts, status)
```

## Simulation

- `simulate_mfgarch`
- `simulate_mfgarch_rv_dependent`
- `simulate_mfgarch_diffusion`

The standard simulator accepts an optional seed, Student-t degrees of freedom,
and correlation between low-frequency and return innovations.

## Low-level upstream kernels

- `sum_tau`
- `sum_tau_fcts`
- `calculate_h_andersen`
- `calculate_p`
- `simulate_r`

## Utilities

- `print_fit_summary`
- `write_simulation_csv`
- `parameter_names`
- `model_parameters`
- `model_to_raw` and `raw_to_model`
