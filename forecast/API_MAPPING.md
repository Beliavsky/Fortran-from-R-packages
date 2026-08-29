# API mapping

This table maps the major numerical `forecast` 9.0.2 workflows to direct Fortran APIs.
R formula/S3/plotting layers are intentionally not mapped one-for-one.

| Upstream R API / kernel | Fortran API | Status |
|---|---|---|
| `BoxCox`, `InvBoxCox` | `boxcox`, `inv_boxcox` | translated |
| `BoxCox.lambda` / Guerrero | `boxcox_lambda_loglik`, `boxcox_lambda_guerrero` | translated |
| `Acf`, `Pacf` | `acf_values`, `pacf_values`, `autocorrelation` | translated |
| `Ccf` | `cross_correlation` | translated numerical layer |
| `fourier`, `seasonaldummy` | `fourier_terms`, `seasonal_dummy` | translated |
| `ma`, differencing helpers | `moving_average`, `difference_series`, `seasonal_difference` | translated |
| `findfrequency` | `findfrequency` | translated heuristic |
| `accuracy` | `accuracy` | translated |
| `dm.test` | `dm_test` | translated; normal tail approximation |
| `meanf` / `mean_model` | `mean_forecast` | translated |
| `rwf`, `naive` | `random_walk_forecast`, `naive_forecast` | translated |
| `snaive` | `seasonal_naive_forecast` | translated |
| `croston`, `croston_model` | `fit_croston`, `forecast_croston` | translated |
| `thetaf`, `theta_model` | `fit_theta`, `forecast_theta` | translated numerical representation |
| `etscalc.c` | `ets_state_forecast`, `ets_update`, `ets_calc` | direct algorithm translation |
| `ets`, `forecast.ets`, `ses`, `holt`, `hw` | `ets_fit`, `ets_auto`, `ets_forecast`, `ets_forecast_simulated`, `ses_fit`, `holt_fit`, `hw_fit` | translated; search, optimized initial states, admissibility/fixed parameters, class-1/2/3 and simulation intervals |
| `simulate.ets` | `ets_simulate` | translated |
| `Arima` | `arima_fit`, `arima_refit` | translated; CSS approximation plus stationary/diffuse Gaussian ML, numeric xreg and missing-observation filtering |
| `auto.arima` | `auto_arima` | upstream-ordered stepwise traversal, xreg, IC/start/model-budget controls, fixed d/D, approximation truncation + full ML refit |
| `forecast.Arima` | `arima_forecast` | translated with ARMA+differencing impulse-response variance |
| `simulate.Arima` | `arima_simulate` | stationary/integrated ARIMA simulation with drift and optional xreg effects |
| `arima.errors`, `arimaorder` | `arima_errors`, `arima_order` | translated |
| `ndiffs` | `ndiffs` | translated using bundled `urca` KPSS |
| `nsdiffs` | `nsdiffs` | translated with MSTL seasonal strength and OCSB option |
| `ocsb.test` | `ocsb_statistic`, `ocsb_critical_value` | translated numerical layer |
| `arfima` | `arfima_fit`, `arfima_forecast` | translated using bundled `fracdiff` |
| `nnetar` | `nnetar_fit`, `nnetar_forecast` | translated using bundled `nnet` |
| `makeBATSMatrices.cpp`, `calcBATS.cpp` | `bats_make_w`, `bats_make_g`, `bats_make_f`, `state_space_filter` | translated |
| `makeTBATSMatrices.cpp`, `calcTBATS.cpp` | `tbats_make_w`, `tbats_make_g`, `tbats_make_f`, `state_space_filter` | translated |
| `bats` / `forecast.bats` | `bats_fit`, `bats_auto`, `bats_refit`, `bats_forecast` | translated fitting/selection/refit path with ARMA-error optimization and forced controls |
| `tbats` / `forecast.tbats` | `tbats_fit`, `tbats_fit_real`, `tbats_auto`, `tbats_refit`, `bats_forecast` | translated fitting/selection/refit, non-integer seasons and harmonic search |
| `dshw` | `dshw_fit`, `dshw_forecast` | translated |
| `mstl`, `seasadj`, `seasonal`, `trendcycle`, `remainder` | `mstl_decompose`, `seasadj`, `seasonal_component`, `trendcycle`, `remainder_component` | self-contained LOESS/STL + iterative MSTL; not bit-for-bit `stats::stl` |
| `tslm` | `tslm_fit` | matrix numerical layer translated |
| `forecast.lm` | `regression_forecast` | translated |
| trend/season formula terms | `trend_season_matrix` | explicit-matrix replacement |
| `CV.lm` | `regression_cv_stats` | translated numerical metrics |
| `modelAR` | `modelar_fit`, `modelar_forecast` | translated linear-model path |
| `spline_model`, `splinef` | `spline_model_fit`, `spline_forecast` | covariance/MLE path translated |
| `na.interp` | `na_interp` | translated linear interpolation layer |
| `tsoutliers`, `tsclean` | `ts_outliers`, `ts_clean` | translated robust numerical analogue |
| `MBB`, `bld.mbb.bootstrap` | `moving_block_bootstrap`, `bld_mbb_bootstrap` | translated |
| `baggedETS` | `bagged_ets_forecast` | translated ensemble path with spread summaries |
| `tsCV` | `ts_cv` | translated via procedure callback |
| `checkresiduals` Ljung-Box computation | `ljung_box` | translated statistic; no plotting/printing |
| `modeldf` | `model_df_arima`, `model_df_ets`, `model_df_bats` | translated |
| `CVar`, general `baggedModel` | `cvar`, `bagged_forecast` | translated numerical callback APIs |
| `stlm`, `stlf` | `stlf_forecast` plus `mstl_decompose`/model callbacks | translated numerical reseasonalization path |
| calendar helpers | `month_days_sequence`, `business_days_sequence`, `easter_gregorian`, `easter_effect` | translated numerical helpers |
| future Fourier/seasonal helpers | `fourier_terms_multi`, seasonal extrapolation helpers | translated |
| tapered correlation helpers | `tapered_acf`, tapered PACF, `linear_process_bootstrap`, `tapered_correlation_ci` | translated including bootstrap confidence intervals |
| `splinef` GCV selection | spline fitting APIs | self-contained GCV path translated |
| all `gg*`, `autoplot`, `plot.*`, `GeomForecast` | — | deliberately skipped plotting |
