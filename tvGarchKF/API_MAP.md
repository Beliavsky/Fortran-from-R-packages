# API map

| Upstream R/C routine | Modern Fortran API | Status |
|---|---|---|
| `polynomial` | `polynomial_values`, `evaluate_tv_function` | Exact numerical translation |
| `NoLineal` | `nonlinear_values`, `evaluate_tv_function` | Exact when coefficient/exponent lengths agree |
| `trigonometric` | `trigonometric_values`, `evaluate_tv_function` | Exact basis/coefficient behavior |
| R `parse`/`eval` argument | `argument_kind` or `custom_argument` callback | Typed replacement |
| `checkInput` | validation in `evaluate_tv_function` | Implemented with status codes |
| C `tvGarch1ab` | `tvgarch_kalman_filter` | Unified exact recursion |
| C `tvGarch1ab2` | `tvgarch_kalman_filter` | Unified exact recursion |
| `tvGarchKalmanLoglike` | `tvgarch_kalman_loglike` | Source objective preserved |
| `tvGarchKalmanFit` | `tvgarch_kalman_fit` | Implemented; Nelder-Mead optimizer |
| `tvGarchKalmanPrint` | `tvgarch_kalman_print`, filter result fields | Implemented without printing/plots |
| `tvGarch_Sim` | `tvgarch_simulate` | Implemented, including supplied innovations |
| `tvParameter` | `tv_parameter` | Implemented using supplied fGarch port |
| `indipsa` dataset | retained in upstream snapshot | Not converted to a compiled Fortran array |
| plotting calls | omitted | Non-computational |

Fortran is case-insensitive, so calls such as `tvGarchKalmanFit` resolve to the
exported `tvgarch_kalman_fit` spelling only when underscores are included;
idiomatic underscore-separated names are used throughout the public API.
