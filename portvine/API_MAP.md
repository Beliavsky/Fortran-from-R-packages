# API coverage map

| R `portvine` operation | Fortran counterpart | Coverage |
|---|---|---|
| `default_garch_spec()` | `make_portvine_spec()` | Typed numerical counterpart |
| `marginal_settings()` | `marginal_settings_type`, `make_marginal_settings()` | Direct settings counterpart |
| `vine_settings()` | `vine_settings_type`, `make_vine_settings()` | Direct settings counterpart |
| `est_var()` | `est_var()` | R type-7 empirical quantile |
| `est_es(method="mean")` | `est_es(..., risk_es_mean)` | Direct |
| `est_es(method="median")` | `est_es(..., risk_es_median)` | Direct |
| `est_es(method="mc")` | `est_es(..., risk_es_mc)` | Direct Monte Carlo integration |
| internal `dvine_ordering()` | `greedy_dvine_order()` | Native partial-correlation implementation |
| `cond_dvine1_cpp()` | `conditional_dvine_sample()` with one value | Rosenblatt-scale implementation |
| `cond_dvine2_cpp()` | `conditional_dvine_sample()` with two values | Raw or conditional-quantile inputs |
| internal `rcondvinecop()` | `conditional_dvine_sample()` | One or two leading D-vine variables |
| internal `estimate_marginal_models()` | `fit_rolling_marginals()` | Rolling matrix-first counterpart |
| internal `estimate_dependence_and_risk()` | `fit_vine_windows()` plus `estimate_risk_roll()` | Integrated native workflow |
| `estimate_risk_roll()` | `estimate_risk_roll()` | Main high-level counterpart |
| `roll_residuals()` | retained arrays in `marginal_window_result` | Direct access without R class extraction |
| `risk_estimates()` | `portvine_roll_result%overall`, `%conditional` | Typed array access |
| `fitted_marginals()` | `portvine_roll_result%marginal` | Typed array access |
| `fitted_vines()` | `%dvine` or `%cvine` | Typed array access |

## Result dimensions

- `overall(measure, alpha, forecast_time)`
- `conditional(measure, alpha, condition_case, forecast_time)`
- `conditional_value(condition_variable, condition_case, forecast_time)`
- `weights(asset, vine_window)`

The final conditional case is the observed residual case. Its
`condition_level` is stored as `-1`; preceding entries correspond to `cond_u`.

## Adapted or unavailable interfaces

- R's `rvine` mode uses arbitrary regular-vine structure selection. The native
  dependency currently provides C-vines and D-vines; `vine_cvine` is the
  supported approximation for the unrestricted mode.
- Marginal ARMA and GARCH parameters are fitted in two numerical stages.
- Asset names, dates, data frames, S4 accessors, and printed summaries are
  replaced by indices and typed arrays.
- Parallel `future.apply` execution is serial.
- Plotting, vignettes, and packaged `.rda` data are not translated.
