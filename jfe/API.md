# API

All procedures are available through:

```fortran
use jfe
```

The working real kind is `dp = kind(1.0d0)`.

## Returns and risk

- `return_annualized(r, scale [, geometric, status])`
- `sharpe_ratio(r [, rf, alpha, risk_method, status])`
- `sharpe_ratio_annualized(r [, rf, alpha], scale [, geometric, risk_method, status])`
- `adjusted_sharpe_ratio(r [, rf], scale [, risk_method, alpha, status])`
- `downside_deviation(r [, mar, method, potential, status])`
- `downside_potential(r [, mar, status])`
- `upside_risk(r [, mar, method, statistic, status])`
- `omega_sharpe_ratio(r [, mar, status])`
- `sortino_ratio(r [, mar, status])`
- `prospect_ratio(r, mar [, status])`
- `volatility_skewness(r [, mar, statistic, status])`
- `m2_sortino(ra, rb [, mar], scale [, status])`
- `value_at_risk(r [, alpha, status])`
- `expected_shortfall(r [, alpha, status])`

Risk-method constants are `risk_stddev`, `risk_var`, and `risk_es`.
Downside methods are `downside_full` and `downside_subset`.

## Drawdown-based measures

- `drawdowns(r [, geometric, status])`
- `max_drawdown(r [, geometric, invert, status])`
- `drawdown_peak(r [, status])`
- `ulcer_index(r [, status])`
- `pain_index(r [, status])`
- `pain_ratio(r [, rf], scale [, status])`
- `calmar_ratio(r, scale [, status])`
- `sterling_ratio(r, scale [, excess, status])`
- `burke_ratio(r [, rf], scale [, modified, status])`
- `martin_ratio(r [, rf], scale [, status])`

`drawdowns`/`max_drawdown` use decimal returns. `drawdown_peak`, Ulcer, Pain,
and the Burke episode calculation preserve the upstream percentage-return
convention.

## Relative performance

- `active_premium(ra, rb, scale [, geometric, status])`
- `tracking_error(ra, rb, scale [, status])`
- `information_ratio(ra, rb, scale [, status])`
- `capm_jensen_alpha(ra, rb [, rf], scale [, status])`
- `appraisal_ratio(ra, rb [, rf], scale [, method, status])`
- `treynor_ratio(ra, rb [, rf], scale [, modified, status])`

Appraisal constants are `appraisal_standard`, `appraisal_modified`, and
`appraisal_alternative`.

## Distributional and summary measures

- `bernardo_ledoit_ratio(r [, status])`
- `d_ratio(r [, status])`
- `kelly_ratio(r [, rf, status])`
- `skewness_kurtosis_ratio(r [, status])`
- `mean_absolute_deviation(r [, status])`
- `skewness(r [, method, status])`
- `kurtosis(r [, method, status])`
- `table_annualized_returns(r, scale [, rf, geometric])`

`table_annualized_returns` returns an `annualized_summary` containing
`annualized_return`, `annualized_sd`, `annualized_sharpe`, and `status`.

## Durbin h

```fortran
result = durbin_h(residuals, n_fitted, n_coefficients, lag_variance)
```

The `durbin_h_result` fields are `statistic`, `p_value`, `durbin_watson`, and
`status`. `lag_variance` is the estimated variance of the lagged dependent
variable coefficient from the fitted model.
