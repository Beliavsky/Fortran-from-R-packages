# API mapping

This table maps computational elements of the R package to the translated
Fortran API. R classes, methods, plotting, labels, and formula parsing are not
mapped.

| R or dependency function | Modern Fortran procedure |
|---|---|
| `gogarch(..., estby="ica")` | `fit_gogarch_ica` or `fit_gogarch(...,"ica")` |
| `gogarch(..., estby="mm")` | `fit_gogarch_mm` or `fit_gogarch(...,"mm")` |
| `gogarch(..., estby="nls")` | `fit_gogarch_nls` or `fit_gogarch(...,"nls")` |
| `gogarch(..., estby="ml")` | `fit_gogarch_ml` or `fit_gogarch(...,"ml")` |
| `goinit` | `initialize_gogarch` |
| `Rd2` | `rd2` |
| `UprodR` | `uprod_r` |
| `Umatch` | `umatch` |
| `unvech` | `unvech` |
| implicit inverse operation | `vech` |
| `cora` | `cora` |
| `gonls` | `gonls_objective` |
| `gotheta` | `gogarch_from_angles` |
| `gollh` | `gogarch_negloglik` |
| `fastICA` dependency | `fastica` |
| `garchFit(~garch(p,q))` | `fit_garchpq` or `fit_univariate` with `model="garch"` |
| `garchFit` APARCH formula | `fit_univariate` with `model="aparch"` |
| `cond.dist="norm"` | `distribution="norm"` |
| `cond.dist="snorm"` | `distribution="snorm"` |
| `cond.dist="std"` | `distribution="std"` |
| `cond.dist="sstd"` | `distribution="sstd"` |
| `cond.dist="ged"` | `distribution="ged"` |
| `cond.dist="sged"` | `distribution="sged"` |
| factor `h.t` values | `garch11_fit%variance` and `gogarch_fit%factor_variance` |
| conditional covariance methods | `gogarch_fit%covariance`, `conditional_variances`, `conditional_correlations` |
| residual methods | `standardized_residuals` |
| prediction methods | `forecast_univariate`, `forecast_gogarch` |
| fitted simulation | `simulate_fitted_gogarch` |
| coefficient methods | `factor_coefficients`, `factor_coefficients_full` |

The retained type name `garch11_fit` is backward-compatible with release
0.1.0, but the type now stores arbitrary ARCH, leverage, and GARCH lag vectors
and therefore represents all supported univariate specifications.
