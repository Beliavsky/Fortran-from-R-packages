# API

All public APIs are re-exported by:

```fortran
use ufrisk
```

Real values use `dp` from `kind_mod`.

## Main types

### `ufrisk_options`

Controls `varcast`.

Important fields:

- `model`: one of `ufrisk_model_sgarch`, `ufrisk_model_lgarch`, `ufrisk_model_egarch`, `ufrisk_model_aparch`, `ufrisk_model_figarch`, or `ufrisk_model_filgarch`
- `distribution`: `ufrisk_distribution_normal` or `ufrisk_distribution_student`
- `smooth`: `ufrisk_smooth_none` or `ufrisk_smooth_lpr`
- `arch_order`, `garch_order`
- `n_out`
- `var_confidence`, `es_confidence`
- `smoothing_order`, normally 1 or 3
- `smoothing_algorithm`, matching the translated smoots algorithms
- `fractional_p_min`, `fractional_p_max`, `fractional_q_min`, `fractional_q_max` for long-memory scale-error selection
- `truncation_lag`
- `max_fit_iterations`, `fit_tolerance`

### `ufrisk_result`

Contains:

- model and distribution identifiers
- estimated in-sample mean
- original and centered in/out returns
- in-sample and forecast volatility
- nonparametric scale estimates
- VaR at the ES level and VaR level
- ES forecasts
- Student-t degrees of freedom
- retained `rugarch_fit_t`, `arma_model`, or `fracdiff_model`
- retained short- or long-memory smoothing result
- status and message

The fields `var_tail_probability` and `es_tail_probability` deliberately store tail probabilities, matching the R package result object after it transforms the input confidence levels.

## Forecasting

```fortran
result = varcast(prices [, options])
```

`prices` must be a positive chronological price vector. The routine computes log returns internally, estimates the model on all but the final `n_out` returns, and produces rolling one-step forecasts for those observations.

## Backtesting

```fortran
traffic = trafftest(result)
coverage = covtest(loss, value_at_risk, tail_probability)
losses = lossfunc(loss, expected_shortfall [, beta])
```

Result types are:

- `ufrisk_traffic_result`
- `ufrisk_coverage_result`
- `ufrisk_loss_result`

## Lower-level translated algorithms

```fortran
coefficients = arfilt_coefficients(ar, ma, d, k)
degrees_freedom = estimate_student_df(standardized_residuals)
call long_memory_smooth(y, result, ...)
```

`arfilt_coefficients` translates the hidden R `arfilt` routine used by both Log-GARCH variants.

## Status values

- `ufrisk_ok`
- `ufrisk_invalid_input`
- `ufrisk_smoothing_failed`
- `ufrisk_model_fit_failed`
- `ufrisk_numerical_failure`
- `ufrisk_no_violations`
