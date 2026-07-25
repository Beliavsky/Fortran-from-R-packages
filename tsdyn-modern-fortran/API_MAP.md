# API map

This map distinguishes implemented numerical translations from partial analogues and exclusions.

## Implemented

| R/package area | Fortran procedures |
|---|---|
| `linear`, `lineVar`, `linear.sim`, `linear.boot` | `fit_ar`, `fit_var`, `simulate_ar`, `simulate_var`, `residual_bootstrap_ar`, `bootstrap_var` |
| `charac_root`, `lags.select` | `ar_roots`, `select_ar_order`, `select_var_lag` |
| `setar`, `setar.sim`, `setar.boot`, `selectSETAR`, regime extraction | `fit_setar`, `simulate_setar`, `residual_bootstrap_ar`-style innovations through simulation, `select_setar_orders`, `setar_regimes` |
| `lstar`, `selectLSTAR`, sigmoid transition | `fit_lstar`, `simulate_lstar`, `forecast_lstar`, `lstar_transition` |
| `llar`, `llar.fitted`, `llar.predict` | `llar_fit_curve`, `llar_fitted`, `llar_predict` |
| `VAR.sim`, `VAR.boot`, prediction | `simulate_var`, `bootstrap_var`, `forecast_var` |
| `VECM`, `VECM.sim`, prediction, `VARrep` | `fit_vecm`, `simulate_vecm`, `forecast_vecm`, `vecm_var_coefficients` |
| `rank.select`, rank statistics | `select_vecm_rank`, `johansen_statistics` |
| `TVAR`, `TVAR.sim`, prediction, regime | `fit_tvar`, `simulate_tvar`, `forecast_tvar`, `tvar_regimes` |
| `TVECM`, `TVECM.sim`, prediction, regime | `fit_tvecm`, `simulate_tvecm`, `forecast_tvecm`, `tvecm_regimes` |
| linear IRF and FEVD | `impulse_response_var`, `impulse_response_vecm`, `fevd_var` |
| nonlinear/regime IRF and GIRF | `regime_irf_setar`, `regime_irf_tvar`, `girf_setar`, `girf_tvar`, `girf_tvecm` |
| `delta`, `delta.lin`, `delta.test`, `delta.lin.test` | `delta_statistic`, `delta_linear_statistic`, `delta_shuffle_test`, `delta_linearity_test` |
| `BBCTest` | `bbc_unit_root_test` |
| `KapShinTest` | `kapshin_test` |
| `resample_vec` | `resample_vector`, `block_resample_matrix` |
| `predict_rolling` | `rolling_forecast_ar`, `rolling_forecast_setar`, `rolling_forecast_var` |
| `accuracy_stat`, `MAPE`, `mse` | `compute_accuracy` and `accuracy_metrics` |

## Partial numerical analogues

- `setarTest` and `TVAR.LRtest`: threshold search and LR statistics are available through `setar_lr_statistic` and `tvar_lr_statistic`; complete package-specific bootstrap distributions are not.
- `TVECM.HStest` and `TVECM.SeoTest`: `tvecm_lr_statistic` supplies a fitted linear-versus-threshold LR calculation and selected threshold, not the complete published bootstrap algorithms.
- `star`: the two-regime LSTAR model is implemented; generic multiple smooth transitions, model mutation, analytical gradients, and all starting-value machinery are not.
- Johansen routines: rank-restricted ML estimation and trace/max-eigen statistics are implemented directly, but exact `urca` normalization and deterministic-case conventions are not promised.
- Local-linear AR: the numerical model is implemented; plotting/data-frame methods and the original C memory/index layout are not.

## Excluded

- `aar`/`mgcv` additive AR models
- `nnetTs` and `selectNNET`
- Full Hansen 1999, Hansen-Seo 2002, and Seo 2006 bootstrap test implementations
- `VECM_symbolic`
- Plotting, GUI, S3 dispatch, formulas, date/time-series metadata, `tidy`, `toLatex`, and print/summary methods
- External-package wrappers and exact R RNG/optimizer compatibility
