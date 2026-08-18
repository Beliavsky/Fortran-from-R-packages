# API map

The upstream `miscTools` 0.6-30 NAMESPACE exports 25 names.

## Fully translated computational utilities

| R export | Fortran mapping |
|---|---|
| `checkNames` | `check_names` |
| `coefTable` | `coef_table` |
| `colMedians` | `col_medians`, `col_medians_3d`, `col_medians_4d` |
| `rowMedians` | `row_medians` |
| `ddnorm` | `ddnorm` |
| `insertCol` | `insert_col` |
| `insertRow` | `insert_row` |
| `isSemidefinite` | `is_semidefinite` |
| `semidefiniteness` | `semidefiniteness` |
| `quasiconcavity` | `quasiconcavity` |
| `quasiconvexity` | `quasiconvexity` |
| `rSquared` | `r_squared` |
| `stdEr` | `std_er` from a supplied covariance matrix |
| `symMatrix` | `sym_matrix` |
| `triang` | `triang` |
| `vecli` | `vecli` |
| `vecli2m` | `vecli2m` |
| `veclipos` | `veclipos` |

## Consolidated/model-interface equivalents

- `nObs`: `n_obs_matrix` returns the number of observations represented by
  a numerical design/data matrix.
- `nParam`: `n_param_vector` returns the number of supplied coefficients.

R's generic/S3 dispatch onto arbitrary fitted model objects is intentionally
not reproduced.

## Intentionally omitted R-only/presentation APIs

- `compPlot` — plotting
- `histDens` — histogram/density plotting
- `summarizeDF` — console/file reporting plus `digest` hashes
- `sumKeepAttr` — specifically preserves arbitrary R attributes on a scalar
- `margEff` — the supplied package has only a generic and an erroring default
  method; there is no computational marginal-effects implementation to port

The ordinary numerical sum corresponding to `sumKeepAttr` is directly
available as Fortran's intrinsic `sum`.
