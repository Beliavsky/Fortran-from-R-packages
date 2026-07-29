# API

All public interfaces are re-exported by module `apt`.

## Model constants

- `apt_linear = 0`
- `apt_tar = 1`
- `apt_mtar = 2`

## Threshold cointegration

### `ci_tar_fit(y, x, result, model, lag, threshold, start_index)`

Fits the long-run regression `y = intercept + slope*x + error`, then estimates

`delta(error_t) = rho_pos*I_t*error_(t-1) + rho_neg*(1-I_t)*error_(t-1) + lagged error differences`.

`model=apt_tar` sets `I_t` from the lagged error. `model=apt_mtar` sets it from the lagged error difference. `start_index` optionally imposes a common estimation start for lag comparisons.

Original-name generic: `ciTarFit`.

### `ci_tar_lag(y, x, result, model, max_lag, threshold, adjust)`

Fits lags `0:max_lag`, records SSE/AIC/BIC and Ljung-Box p-values, and selects the AIC- and BIC-minimizing lags. With `adjust=.true.`, all candidate regressions use a common sample.

Original-name generic: `ciTarLag`.

### `ci_tar_threshold(y, x, result, model, lag, trim_fraction)`

Sorts candidate TAR or MTAR threshold variables, removes the requested fraction from each tail, fits every remaining threshold, and selects the first SSE minimum.

Original-name generic: `ciTarThd`.

## Error-correction models

### `ecm_symmetric_fit(y, x, result, lag)`

Fits two ECM equations for `delta(x)` and `delta(y)` with an intercept, lagged changes in both variables, and one lagged error-correction term.

Original-name generic: `ecmSymFit`.

### `ecm_asymmetric_fit(y, x, result, lag, split, model, threshold)`

Fits two ECM equations with positive and negative error-correction terms. With `split=.true.`, each lagged price change is also split into positive and negative parts.

Original-name generic: `ecmAsyFit`.

### `ecm_asymmetry_tests(model, result)`

Returns paired equation tests for:

1. equilibrium-adjustment path symmetry;
2. x-lag Granger exclusion;
3. y-lag Granger exclusion;
4. per-lag positive/negative symmetry when split;
5. cumulative positive/negative symmetry when split.

Original-name generic: `ecmAsyTest`.

### `ecm_diagnostics(model, result)`

Returns R-squared, adjusted R-squared, equation F statistic, Durbin-Watson statistic and approximate p-value, AIC, BIC, and Ljung-Box p-values at lags 4, 8, and 12.

Original-name generic: `ecmDiag`.

## Statistical foundation

- `fit_ols`
- `linear_f_test`
- `zero_coefficient_f_test`
- `ljung_box_test`
- `durbin_watson_test`
- `normal_cdf`
- `student_t_cdf`
- `f_cdf`
- `chi_square_cdf`

## Result types

- `regression_result`
- `hypothesis_result`
- `ci_tar_fit_result`
- `ci_tar_lag_result`
- `ci_tar_threshold_result`
- `ecm_fit_result`
- `ecm_diagnostics_result`
- `ecm_asymmetry_test_result`

Every high-level result contains an integer `status`; zero means success.
