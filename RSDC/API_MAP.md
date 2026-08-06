# RSDC to Fortran API map

| RSDC routine | Fortran routine/type | Status |
|---|---|---|
| `rsdc_hamilton` | `rsdc_hamilton` / `rsdc_filter_result` | Formula-equivalent filter and smoother |
| `rsdc_likelihood` | `rsdc_negative_log_likelihood` | Implemented for const, noX, TVTP |
| `rsdc_estimate` | `rsdc_estimate` / `rsdc_model` | Implemented; jDE + pattern search |
| `rsdc_starts` | `rsdc_make_starts` / `rsdc_starts_result` | Implemented |
| `rsdc_forecast` | `rsdc_forecast_path` / `rsdc_forecast_result` | Implemented |
| `rsdc_forecast_ahead` | `rsdc_forecast_ahead` | Implemented |
| `rsdc_simulate` | `rsdc_simulate`; `rsdc_simulate_fixed` | TVTP and fixed-P simulation |
| `rsdc_viterbi` | `rsdc_viterbi_path` | Implemented |
| `rsdc_minvar` | `rsdc_minvar` / `rsdc_portfolio_result` | Implemented |
| `rsdc_maxdiv` | `rsdc_maxdiv` / `rsdc_portfolio_result` | Implemented |
| `rsdc_bootstrap` | `rsdc_parametric_bootstrap` / `rsdc_bootstrap_result` | Draws, covariance, SEs, and percentile intervals |
| `rsdc_corr_bands` | `rsdc_corr_bands` / `rsdc_bands_result` | Implemented |
| score/OPG/sandwich helpers | `rsdc_scores`, `rsdc_robust_vcov` | Implemented |
| transition diagnostics | `rsdc_diagnostics` | Implemented |
| canonical partial-correlation helpers | `partial_to_correlation`, `correlation_to_partial` | Implemented |
| `coef`, `logLik`, `vcov`, `confint`, print/summary | fields of `rsdc_model` | R generic wrappers omitted |
| broom and plotting methods | none | Omitted as non-computational UI |
| R parallel backend | none | Omitted; caller may parallelize independent fits |

## Principal result fields

`rsdc_model` contains natural packed parameters, TVTP coefficients, regime
correlations, representative transition matrix, covariance slices, in-sample
and optional OOS log likelihoods, convergence status, and optional covariance
and standard-error arrays.
