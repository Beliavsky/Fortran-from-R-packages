# Source coverage

## Translated computational files

- `R/mod_duration.R`
- `R/convexity.R`
- `R/total_return.R`
- the `na.locf0` and percent-conversion steps in `R/get_yields.R`

## Replaced with native Fortran representations

- R vectors and `xts` columns: rank-1 arrays
- multi-column `xts` objects: rank-2 `(time, series)` arrays
- R `NA_real_`: IEEE quiet NaN
- R exceptions: optional integer status and diagnostic message

## Omitted files or branches

- FRED download in `R/get_yields.R`
- `R/tibble_to_xts.R`
- `R/xts_to_tibble.R`
- `R/dependencies.R`
- plotting and package/vignette infrastructure

All upstream source files remain in `provenance/upstream-source`.
