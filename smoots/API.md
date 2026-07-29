# API reference

All public procedures are available from:

```fortran
use smoots
```

## Smoothing and estimation

### `fixed_gsmooth(y, v, p, mu, b, bb, result)`

Fixed-bandwidth local polynomial estimator corresponding to R `gsmooth`.

- `v`: derivative order
- `p`: polynomial order; `p > v` and `p-v` must be odd
- `mu`: kernel smoothness, 0 through 3
- `b`: bandwidth fraction, between 0 and 0.5
- `bb`: 0 for fixed boundary windows, 1 for boundary-adaptive windows

`result%weights` contains the complete boundary/interior equivalent-kernel weight system.

### `fixed_knsmooth(y, mu, b, bb, result)`

Fixed-bandwidth Nadaraya-Watson kernel smoother corresponding to R `knsmooth`.

### `tsmooth(y, p, mu, cf_method, inflation, b_start, enlarged_variance, bb, boundary_cut, method, result)`

Fully configurable iterative plug-in trend estimator.

Variance-factor constants:

- `sm_cf_np`
- `sm_cf_ar`
- `sm_cf_ma`
- `sm_cf_arma`

Inflation constants:

- `sm_infl_opt`
- `sm_infl_naive`
- `sm_infl_variance`

Smoothing methods:

- `sm_method_lpr`
- `sm_method_kernel`

### `msmooth(y, result [, p, mu, b_start, algorithm, method])`

Convenience interface implementing all original named algorithms:

| Algorithm | Variance estimate | Inflation | Enlarged variance bandwidth |
|---|---|---|---|
| `A` | nonparametric | optimal | yes |
| `B` | nonparametric | naive | yes |
| `O` | nonparametric | optimal | no |
| `N` | nonparametric | naive | no |
| `OA` | AR | optimal | no |
| `NA` | AR | naive | no |
| `OM` | MA | optimal | no |
| `NM` | MA | naive | no |
| `OAM` | ARMA | optimal | no |
| `NAM` | ARMA | naive | no |

### `dsmooth(y, result [, d, mu, pilot_p, pilot_b_start, b_start])`

Data-driven derivative estimation. `d` is 1 or 2.

### `conf_bounds(object [, confidence, result, parametric_order])`

Returns asymptotically unbiased estimates, pointwise normal bounds, and a parametric comparison fit.

### Low-level smoothers

```fortran
call gsmooth(y, v, p, mu, bandwidth, bb, estimate, weights, status)
call knsmooth(y, mu, bandwidth, bb, estimate, status)
call local_polynomial_smooth(...)
call kernel_smooth(...)
```

## Long-run variance

```fortran
call lag_window_variance(x, cf0, l0_opt, lg_opt, status)
call estimate_cf0_ar(x, cf0, model, status)
call estimate_cf0_ma(x, cf0, model, status)
call estimate_cf0_arma(x, cf0, model, status)
```

The first procedure implements the original Rcpp Bühlmann recursion.

## ARMA utilities

```fortran
call fit_arma(x, p, q, include_mean, model, status)
call arma_residuals(x, ar, ma, mean_x, residuals, fitted)
call arma_point_forecast(x, residuals, ar, ma, mean_x, h, forecast, status)
call simulate_arma(ar, ma, mean_x, n, burn, series, status, ...)
call ma_infinity(ar, ma, m, coefficients, status)
call information_criterion_matrix(x, pmax, qmax, include_mean, use_bic, matrix, status)
call optimal_order(matrix, p, q, minimize, mask, status)
```

`optimal_order` accepts an optional logical matrix `mask`, which is the typed Fortran replacement for the arbitrary R expression accepted by `optOrd`.

## Forecasting

```fortran
call trend_forecast(model, h, mode, forecast, status)
call normal_forecast(x, p, q, include_mean, h, confidence, result)
call bootstrap_forecast(x, p, q, include_mean, h, simulations, burn, confidence, result [, seed])
call model_forecast(smooth_model, p, q, h, use_bootstrap, confidence, result, ...)
call rolling_backtest(y, k, p, q, use_bootstrap, confidence, result, ...)
```

Trend modes are `sm_trend_linear` and `sm_trend_constant`.

The bootstrap implementation is serial and deterministic when an explicit 64-bit seed is supplied. It retains the original forward-bootstrap sequence: residual resampling, ARMA simulation, parameter re-estimation, proxy residual calculation on the original sample, and comparison of simulated true values with bootstrap forecasts.

## Rescaling

```fortran
rescaled = rescale_derivative(values, first_x, second_x, last_x, order)
```
