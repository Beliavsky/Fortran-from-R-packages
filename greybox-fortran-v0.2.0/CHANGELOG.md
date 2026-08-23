# Changelog

## 0.2.0

- Added two-block beta regression to `alm_fit`, including beta prediction and covariance blocks.
- Added LASSO, RIDGE, ROLE and QUALE ALM objectives.
- Added `alm_fit_occurrence` / `alm_occurrence_model`.
- Added conditional ARIMA-error `alm_fit_arima_errors` / `alm_dynamic_model`.
- Added beta-aware `calm_fit` / `calm_predict`.
- Added `lm_dynamic_fit` with point IC weights and LOWESS-style robust smoothing.
- Added `rmcb_test` with Tukey, normal and general ALM branches.
- Added `dsr_bootstrap`, including intermittent-demand reconstruction.
- Added `aid_fit` and `aid_cat` demand classification.
- Added `test_parity_v02` and kept all v0.1 tests green.

## 0.1.0

- Initial modern Fortran/FPM translation of the core greybox distributions, forecast measures, utilities, association diagnostics, ALM regression, point ICs, stepwise selection, static model averaging, rolling-origin evaluation and recursive regression.
