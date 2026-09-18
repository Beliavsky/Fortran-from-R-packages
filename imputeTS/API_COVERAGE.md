# API coverage

## Basis

Coverage counts distinct exported computational R functions. Plotting, presentation, data, pipe/operator, and internal helper functions are excluded. The 10 deprecated dot-name imputation wrappers are counted because they remain exported computational entry points, but each maps to the same implementation as its underscore-name successor.

**Package status:** substantial  
**Mapped computational functions:** 21 of 21 (100%)

The percentage measures mapped exported computational functions, not exact R interface compatibility.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `na_interpolation` | `imputets_interpolation` | `na_interpolation` | substantial |
| `na_kalman` | `imputets_kalman` | `na_kalman` | partial |
| `na_locf` | `imputets_basic` | `na_locf` | substantial |
| `na_ma` | `imputets_basic` | `na_ma` | substantial |
| `na_mean` | `imputets_basic` | `na_mean` | substantial |
| `na_random` | `imputets_basic` | `na_random` | substantial |
| `na_remove` | `imputets_basic` | `na_remove` | complete |
| `na_replace` | `imputets_basic` | `na_replace` | substantial |
| `na_seadec` | `imputets_seasonal` | `na_seadec` | substantial |
| `na_seasplit` | `imputets_seasonal` | `na_seasplit` | substantial |
| `statsNA` | `imputets_stats` | `stats_na` | substantial |
| `na.interpolation` | `imputets_interpolation` | `na_interpolation` | substantial |
| `na.kalman` | `imputets_kalman` | `na_kalman` | partial |
| `na.locf` | `imputets_basic` | `na_locf` | substantial |
| `na.ma` | `imputets_basic` | `na_ma` | substantial |
| `na.mean` | `imputets_basic` | `na_mean` | substantial |
| `na.random` | `imputets_basic` | `na_random` | substantial |
| `na.remove` | `imputets_basic` | `na_remove` | complete |
| `na.replace` | `imputets_basic` | `na_replace` | substantial |
| `na.seadec` | `imputets_seasonal` | `na_seadec` | substantial |
| `na.seasplit` | `imputets_seasonal` | `na_seasplit` | substantial |

## Excluded exports

The `%>%` operator and all `ggplot_*` / `plotNA.*` exports are R-specific interface or presentation code. Dataset objects are also outside the computational function denominator. Internal `locf`, `ma`, and `apply_base_algorithm` helpers are not exported R functions and are not counted separately.

## Main differences

- Numeric vectors are the public Fortran data model; R `ts`, tibble, data-frame, S3/class attributes, warnings, and columnwise `tryCatch` behavior are not reproduced.
- R `NA`/`NaN` is represented by an IEEE quiet NaN.
- `maxgap < 0` or an absent `maxgap` means unlimited imputation, corresponding to the R default `Inf`.
- `na_interpolation(..., option="spline")` uses a natural cubic spline and does not expose the full `stats::spline` `...` surface.
- `na_kalman` is intentionally partial: the default structural path is a local-linear-trend Kalman filter/smoother; the ARIMA path reuses `forecast` but constrains differencing to zero so fitted values remain aligned with the original sample. Arbitrary R state-space model lists and exact `KalmanSmooth`/`KalmanRun` initialization are not translated.
- Seasonal algorithms reuse `forecast::findfrequency`/STL equivalents from the sibling `forecast` translation. Base-algorithm `...` forwarding is normalized to each translated algorithm's defaults.
- Random imputation uses the Fortran processor RNG, not R's RNG stream.
- `stats_na` returns computational statistics in `na_stats_result`; printing and the bin presentation table are omitted.

## Special-value and fallback behavior

- Seasonal fallback paths preserve the upstream early return: when decomposition/splitting cannot be used, the outer `maxgap` is not reapplied after the selected base algorithm.
- NaNs are the Fortran missing-value marker. Infinite values are ordinary numeric values; this differs from the upstream Rcpp LOCF helper, which tests `R_finite` internally.
