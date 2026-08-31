# API coverage

Upstream: R package `mitools` 2.4.

| Upstream API | Fortran coverage | Notes |
|---|---|---|
| `MIcombine.default` | Implemented | `mi_combine` accepts scalar vectors or parameter-by-imputation arrays plus corresponding variances/covariances. Rubin point estimates, within/between variance combination, finite `df.complete` adjustment, total covariance, and missing-information fraction follow the upstream formulas. |
| `MIcombine.imputationResultList` | Numerically covered | R method dispatch through `coef()`/`vcov()` is omitted. Callers pass extracted numeric estimates and covariance matrices directly. |
| `vcov.MIresult` | Covered by data field | Use `result%variance`. |
| `summary.MIresult` | Implemented numerically | `mi_summary` returns estimate, standard error, Student-t interval, and missing-information fraction. `log_effect=.true.` reproduces the upstream exponentiation and delta-method SE transformation. Printing is left to the caller. |
| `print.MIresult` | Omitted | Formatting-only R interface. |
| `imputationList.default` | Implemented for numeric data | `imputation_list_from_array` stores a numeric `(row,column,imputation)` cube. Heterogeneous R data frames/factors are outside Fortran's numeric kernel. |
| `dim.imputationList` | Implemented | `imputation_dimensions`. |
| `dimnames.imputationList` | Omitted | R naming/interface metadata. |
| `rbind.imputationList` | Implemented for numeric data | `imputation_rbind`; all imputations are bound consistently. |
| `cbind.imputationList` | Implemented for numeric data | `imputation_cbind`; all imputations are bound consistently. |
| `with.imputationList` | Omitted | Arbitrary R expression/function evaluation. Numerical analyses are called explicitly by Fortran clients. |
| `MIextract` | Omitted | R expression/function dispatch over result objects. |
| `update.imputationList` | Omitted | R expression evaluation and dynamic data-frame column creation. |
| `withPV.default` | Computational selection covered | `pv_select` selects one replicate for all plausible-value variables and `pv_materialize` appends those values to a numeric base matrix. Formula rewriting and arbitrary R action evaluation are intentionally omitted. |
| Database-backed imputation lists | Omitted | DBI/RODBC connection management, SQL query generation, and deferred R expressions are external/R-specific interfaces. |
| `updatesInfilter`, `updatesOutfilter`, `getvars` | Omitted | Helpers for database query planning and R expression evaluation rather than standalone numerical algorithms. |

## Parity assessment

The substantive statistical computation in `mitools` is `MIcombine`; it is
translated directly. The remaining differences are primarily dynamic R object,
formula, database, and evaluation infrastructure, which is intentionally outside
the scope requested for computational translation.
