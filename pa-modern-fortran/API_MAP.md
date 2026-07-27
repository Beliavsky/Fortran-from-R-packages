# API map

| Original R routine or method | Modern Fortran procedure | Status |
|---|---|---|
| `brinson()` single period | `fit_brinson_period` | Implemented and tested |
| `brinson()` multiple periods | `fit_brinson_multi` | Implemented and tested |
| `returns(brinson)` | `brinson_period_result%category_effect`, `%aggregate` | Implemented and tested |
| `returns(brinsonMulti, type="arithmetic")` | `summarize_brinson_multi(..., "arithmetic", ...)` | Implemented and tested |
| `returns(brinsonMulti, type="geometric")` | `summarize_brinson_multi(..., "geometric", ...)` | Implemented and tested |
| `returns(brinsonMulti, type="linking")` | `summarize_brinson_multi(..., "linking", ...)` | Implemented and tested |
| `regress()` single period | `fit_regression_period` | Implemented and tested |
| `regress()` multiple periods | `fit_regression_multi` | Implemented and tested |
| `returns(regression)` | `summarize_regression_period` | Implemented and tested |
| `returns(regressionMulti, type=...)` | `summarize_regression_multi` | All three methods implemented and tested |
| R `model.matrix(... - 1)` path | `build_design_matrix` | Implemented for integer categorical and numeric arrays |
| `exposure()` categorical | `categorical_exposure` | Implemented and tested |
| `exposure()` continuous | `continuous_exposure` | Implemented and tested |
| Multi-period categorical exposure | `categorical_exposure_multi` | Implemented and tested |
| Multi-period continuous exposure | `continuous_exposure_multi` | Implemented and tested |
| `.cat.ret` | Internal calculations in `fit_brinson_period` | Implemented |
| `.aggregate`, `.combine` | Attribution summary result types | Implemented numerically |
| S4 classes and accessors | Plain derived types in `pa_types` | R infrastructure replaced |
| Plot methods | None | Excluded |
| `show()` and formatted `summary()` | None | Excluded; numerical content is returned directly |
