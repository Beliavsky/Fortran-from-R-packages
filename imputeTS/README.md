# imputeTS

Modern free-form Fortran translation of the computational core of the R package **imputeTS 3.4** by Steffen Moritz, Sebastian Gatscha, and contributors. Upstream imputeTS is licensed under GPL-3.

The translation targets placement as a top-level sibling directory in `Beliavsky/Fortran-from-R-packages`. It deliberately skips plotting, interactive/presentation behavior, R classes, and datasets.

## Build

From the repository root layout, keep these translated dependencies as siblings:

```text
Fortran-from-R-packages/
  forecast/
  stinepack/
  imputeTS/
```

Then:

```text
cd imputeTS
fpm build
fpm test
```

`fpm.toml` uses sibling path dependencies:

```toml
forecast-fortran = { path = "../forecast" }
stinepack = { path = "../stinepack" }
```

No BLAS, LAPACK, ARPACK, Rcpp, or translated dependency source is copied into this package. The `forecast` sibling supplies the shared frequency/STL/ARIMA machinery and `stinepack` supplies Stineman interpolation.

## Public API

Use module `imputets_api`. It re-exports `dp`, `na_stats_result`, and:

```fortran
na_interpolation(x, option, maxgap)
na_kalman(x, model, smooth, period, maxgap)
na_locf(x, option, na_remaining, maxgap)
na_ma(x, k, weighting, maxgap)
na_mean(x, option, maxgap)
na_random(x, lower_bound, upper_bound, maxgap)
na_remove(x)
na_replace(x, fill, maxgap)
na_seadec(x, algorithm, period, find_frequency, maxgap)
na_seasplit(x, algorithm, period, find_frequency, maxgap)
stats_na(x, bins)
```

Missing values are IEEE quiet NaNs. For the integer `maxgap`, an absent or negative value means unlimited imputation.

## Translation coverage

**Status: substantial -- 21 of 21 (100%)** exported computational R functions have meaningful Fortran mappings. The fraction measures mapped computational functions, not complete R compatibility. Deprecated dot-name wrappers are counted as R entry points but map to the same Fortran procedures as their underscore-name replacements.

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

See `API_COVERAGE.md` for counting details and material compatibility differences.

## Numerical scope

The direct algorithms translate the upstream numerical behavior for global-statistic replacement, LOCF/NOCB, moving-average imputation, uniform random imputation, removal/replacement, `maxgap`, and missing-gap statistics. Stineman interpolation delegates to the translated sibling `stinepack`. Seasonal frequency detection and STL delegate to the translated sibling `forecast`. The Kalman ARIMA branch also uses the sibling `forecast`.

`na_kalman` remains the largest deliberate parity gap. The translated structural mode uses a local-linear-trend state-space model and RTS smoothing instead of reproducing every `stats::StructTS` model variant. The `auto.arima` mode performs shared ARMA model selection with zero differencing for index-aligned imputation and optionally combines forward/reverse filtered values as a two-sided estimate.

For seasonal fallback cases where no usable frequency is available or fewer than two periods exist, the translation follows upstream behavior and invokes the base algorithm directly; the outer `maxgap` argument is therefore not reapplied in that early-return path. The maintained Fortran code treats IEEE NaN as the missing-value marker. Unlike the upstream Rcpp LOCF kernel, infinities are not treated as non-finite missing candidates.

## Tests

`test/test_imputets.f90` deterministically exercises each algorithm family, including the three moving-average weights, all global mean variants, Stineman/spline/linear interpolation, the two Kalman branches, seasonal split/decomposition, random bounds, `statsNA`-style gap summaries, and `maxgap`.

See `VALIDATION.md` for the exact validation environment and limitations.

## Upstream provenance

Retained computational R/C++ source excerpts are under `upstream/` for auditability. Deprecated plotting wrappers are intentionally omitted from the retained `deprecated_defunct.R` excerpt because plotting is outside scope. See `NOTICE.md` and `PROVENANCE.md`.
