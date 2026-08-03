# Translation coverage

## Exported R functions

| R export | Fortran status |
|---|---|
| `fit_volume` | Implemented as `fit_volume` |
| `decompose_volume` | Implemented as `decompose_volume` |
| `forecast_volume` | Implemented as `forecast_volume` |
| `generate_plots` | Omitted: plotting-only |

## Internal computational routines

| R routine/family | Fortran equivalent |
|---|---|
| `spec_volume_model` | `volume_model_spec`, `initialize_volume_spec` |
| `specify_uniss` | default/specification setup in `fit_volume` |
| `uniss_kalman` | `uniss_kalman` |
| `uniss_em_alg` | ordinary branch of `fit_volume` |
| `uniss_em_alg_acc` | accelerated branch of `fit_volume` |
| `smooth_volume_model` | analysis mode of `decompose_volume` |
| `forecast_volume_model` | forecast mode of `decompose_volume` |
| `clean_data` | `clean_volume_data` for matrix input |
| MAE/MAPE/RMSE helpers | `compute_error_metrics` |
| parameter-list validation | typed masks plus `parameters_valid` |

`intraday_xts_to_matrix` is not compiled because Fortran has no built-in `xts` or
calendar-index object. Callers supply the rectangular matrix directly.
