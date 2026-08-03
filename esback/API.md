# API

All real inputs use `real(dp)` from module `esback`.

## Backtests

### `er_backtest(r, q, e, result [, s, b])`

McNeil-Frey exceedance-residual bootstrap test. `b` defaults to 1000. When `s` is supplied, standardized residual p-values are also calculated.

### `cc_backtest(r, q, e, alpha, result [, s, hommel])`

Nolde-Ziegel conditional-calibration test. Without `s`, returns the simple two-dimensional test. With `s`, also returns the general scalar/two-sided and four-instrument/one-sided tests.

### `esr_backtest(r, e, alpha, version, result [, q, b, options])`

Bayer-Dimitriadis ESR test:

- `version=1`: strict ESR, `r ~ e | e`
- `version=2`: auxiliary ESR, `r ~ q | e`; `q` is required
- `version=3`: strict-intercept ESR, `(r-e) ~ e | 1`

Set `b > 0` for iid bootstrap p-values.

## Joint VaR-ES regression

### `esreg_fit(xq, xe, y, alpha, fit [, options, compute_covariance])`

Fits linear conditional VaR and ES equations using the Fissler-Ziegel loss with the esback configuration `g1=2`, `g2=1`.

### `esreg_covariance(...)`

Computes the asymptotic sandwich covariance matrix.

### `esr_loss(r, q, e, alpha [, g1, g2, return_mean])`

Evaluates the joint VaR/ES loss.

### Nuisance estimators

- `quantile_regression`
- `density_quantile_function`
- `conditional_mean_sigma`
- `conditional_truncated_variance`
- `cdf_at_quantile`

## Options

`type(esreg_options)` controls:

- sparsity method: `sparsity_iid`, `sparsity_nid`
- truncated variance: `sigma_ind`, `sigma_scl_n`, `sigma_scl_sp`
- bandwidth: Bofinger, Chamberlain, Hall-Sheather
- misspecification adjustment
- optimization limits, tolerance, starts, and seed

## Results

- `er_backtest_result`
- `cc_backtest_result`
- `esr_backtest_result`
- `esreg_fit_result`

Each contains a status code. ESR results retain the fitted coefficients, predictions, covariance matrix, and loss.
