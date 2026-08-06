# API map

| Upstream R/C API | Fortran API | Status |
|---|---|---|
| `fkf()` / `FKF` / `cfkf` | `fkf`, `kalman_filter` | Ported |
| `fks()` / `FKS` / `cfks` | `fks`, `kalman_smooth` | Ported |
| `plot.fkf()` | none | Omitted: plotting |
| `plot.fks()` | none | Omitted: plotting |
| S3 result objects | `type(fkf_result)`, `type(fks_result)` | Replaced with typed results |
| R input checks | `validate_model` and status codes | Ported |
| BLAS/LAPACK SPD inversion | self-contained Cholesky inverse/log determinant | Equivalent |
| R `NA` support | IEEE quiet NaN support | Ported |

## Result-name correspondence

- `att` -> `fkf_result%att`
- `at` -> `fkf_result%at`
- `Ptt` -> `fkf_result%ptt`
- `Pt` -> `fkf_result%pt`
- `vt` -> `fkf_result%vt`
- `Ft` -> `fkf_result%ft`
- `Ftinv` -> `fkf_result%ftinv`
- `Kt` -> `fkf_result%kt`
- `logLik` -> `fkf_result%log_likelihood`
- smoother `ahatt` -> `fks_result%ahatt`
- smoother `Vt` -> `fks_result%vt`
