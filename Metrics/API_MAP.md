# API map

All upstream exported computational functions are represented.

| Upstream R function | Fortran procedure | Notes |
|---|---|---|
| `bias` | `bias` | Vector aggregate |
| `percent_bias` | `percent_bias` | Preserves Inf/NaN zero-denominator behavior |
| `se`, `ae`, `ape`, `sle`, `ll` | same names | Elemental scalar procedures; array calls work automatically |
| `sse`, `mse`, `rmse` | same names | Vector aggregates |
| `mae`, `mdae`, `mape`, `smape` | same names | Vector aggregates |
| `msle`, `rmsle` | same names | Vector aggregates |
| `rae`, `rse`, `rrse` | same names | Relative error measures |
| `explained_variation` | same name | `1 - rse` |
| `mase` | `mase` | Optional seasonal `step_size` |
| `auc` | `auc` | Average tied ranks and Mann-Whitney statistic |
| `logLoss` | `logloss` | Fortran is case-insensitive |
| `precision`, `recall`, `fbeta_score` | same names | Binary integer labels |
| `ce`, `accuracy` | same names | Generic integer, real, and character overloads |
| `ScoreQuadraticWeightedKappa` | `scorequadraticweightedkappa` and `score_quadratic_weighted_kappa` | Direct and idiomatic spellings |
| `MeanQuadraticWeightedKappa` | `meanquadraticweightedkappa` and `mean_quadratic_weighted_kappa` | Direct and idiomatic spellings |
| `f1`, `apk` | same names | Generic integer, real, and character overloads |
| `mapk` | `mapk` | Arrays of `integer_vector`, `real_vector`, or `string_vector` |

Most aggregate routines accept an optional integer `stat` argument. Status
constants are exported by `metrics`:

- `metrics_success`
- `metrics_invalid_size`
- `metrics_invalid_argument`
- `metrics_undefined`
