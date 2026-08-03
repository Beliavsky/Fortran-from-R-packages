# Translation coverage

| R export | Fortran API | Status |
|---|---|---|
| `fit_AR1_Gaussian` | `fit_ar1_gaussian` | translated |
| `impute_AR1_Gaussian` | `impute_ar1_gaussian` | translated |
| `impute_rolling_AR1_Gaussian` | `impute_rolling_ar1_gaussian` | translated |
| `fit_AR1_t` | `fit_ar1_t` | translated |
| `impute_AR1_t` | `impute_ar1_t` | translated |
| `fit_VAR_t` | `fit_var_t` | translated with documented missing-data adaptation |
| `impute_OHLC` | `impute_ohlc` | translated |
| `impute_Vol` | `impute_vol` | translated |
| `plot_imputed` | none | plotting omitted |

Internal computational support includes missing-block discovery, conditional
Gaussian moments, consecutive-pair initializers, outlier detection, Student-t
CDF evaluation, gamma/normal sampling, dense linear solves, Cholesky sampling,
and multivariate-t weighted regression.
