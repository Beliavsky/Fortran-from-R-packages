# API

All real values use `dp = kind(1.0d0)` from `intraday_kinds`.

## Main types

### `type(volume_parameters)`

Contains:

- `a_eta`: daily-component persistence across day boundaries.
- `a_mu`: dynamic-component AR(1) coefficient at every bin.
- `var_eta`: daily innovation variance at day boundaries.
- `var_mu`: dynamic innovation variance at every transition.
- `r`: observation variance.
- `x0(2)`, `v0(2,2)`: initial state mean and covariance.
- `phi(:)`: centered log-seasonal profile, one value per bin.

### `type(parameter_mask)`

Logical selectors for `a_eta`, `a_mu`, `var_eta`, `var_mu`, `r`, `phi`, `x0`, and
`v0`.

### `type(volume_model_spec)`

Contains `fixed`, `initial`, `is_fixed`, and `has_initial`. Call
`initialize_volume_spec(spec, n_bin)` before setting masks and values.

### `type(volume_fit_control)`

- `acceleration`: use the package's accelerated EM scheme; default `.true.`.
- `maxit`: maximum major iterations; default 3000.
- `abstol`: Euclidean parameter-change tolerance; default `1e-4`.
- `save_history`: retain parameter history.
- `verbose`: progress-print level.

### `type(volume_model)`

Returns fitted parameters, fixed-parameter mask, convergence status, iteration count,
final change, message, and optional parameter history.

### `type(volume_decomposition)`

Returns flattened original and fitted signals, daily/dynamic/seasonal/residual
components, MAE/MAPE/RMSE, mode, status, and message.

## Main procedures

```fortran
subroutine fit_volume(data, model, spec, control)
  real(dp), intent(in) :: data(:, :)
  type(volume_model), intent(out) :: model
  type(volume_model_spec), intent(in), optional :: spec
  type(volume_fit_control), intent(in), optional :: control
end subroutine
```

```fortran
subroutine decompose_volume(purpose, model, data, result, burn_in_days)
  character(len=*), intent(in) :: purpose
  type(volume_model), intent(in) :: model
  real(dp), intent(in) :: data(:, :)
  type(volume_decomposition), intent(out) :: result
  integer, intent(in), optional :: burn_in_days
end subroutine
```

`purpose` is `"analysis"`, `"smooth"`, or `"forecast"`. `burn_in_days` applies only
to forecasting.

```fortran
subroutine forecast_volume(model, data, result, burn_in_days)
```

A wrapper around `decompose_volume("forecast", ...)`.

## Lower-level procedures

```fortran
subroutine uniss_kalman(log_data, par, output, smooth)
subroutine uniss_em_update(log_data, current, fixed, updated, status, message)
subroutine clean_volume_data(data, cleaned, kept_columns, status, message)
subroutine simulate_intraday_volume(par, n_day, volume, states, seed, status)
```

Status constants are `intraday_ok`, `intraday_invalid_input`,
`intraday_not_converged`, and `intraday_numerical_failure`.
