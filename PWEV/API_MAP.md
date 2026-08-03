# API map

| R functionality | Fortran counterpart | Notes |
|---|---|---|
| `PWEV(Data, SplitR)` | `pwev_fit(data, split_ratio, result, status, control)` | Full high-level workflow. |
| `WeightedEnsemble(df, Method="PSO")` | `pso_ensemble_weights(actual, forecasts, control, result)` | Implements the PSO path used by PWEV. |
| Weighted fitted result | `matmul(forecasts, weights)` | Weights remain in `[0,1]` and are not normalized. |
| Accuracy helper | `pwev_metric_vector`, `pwev_accuracy_table` | RMSE, MAPE, MAE, RRSE, MDAE, RMSLE, RAE, SMAPE, squared correlation. |
| `ugarchfit` sGARCH | `rugarch` `fit_garch11` through `fit_pwev_base_models` | Attached dependency. |
| `ugarchfit` gjrGARCH | `rugarch` `fit_gjrgarch11` through `fit_pwev_base_models` | Attached dependency. |
| `ugarchfit` iGARCH | `rugarch` `fit_igarch11` through `fit_pwev_base_models` | Attached dependency. |
| `umemfit("MEM", "NO", ...)` | `rumidas` `umemfit` through `fit_pwev_base_models` | Attached dependency. |
| R list result | `type(pwev_result)` | Typed arrays and status fields. |

The public Fortran procedure is named `pwev_fit`. Fortran is case-insensitive, so a procedure named `PWEV` cannot coexist cleanly with the top-level module named `pwev`.
