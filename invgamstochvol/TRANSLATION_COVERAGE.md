# Translation coverage

Upstream package: `invgamstochvol` 1.0.0.

| R export | Fortran API | Status |
|---|---|---|
| `ourgeo` | `ourgeo` | Complete |
| `lik_clo` | `lik_clo` | Complete |
| `DrawK0` | `draw_k0` | Complete |

## Internal computational routines

| Upstream routine | Fortran equivalent | Status |
|---|---|---|
| `lrfact` | `log_rising_factorial` | Complete |
| `CalcuLogfac` | `build_log_factorials` | Complete |
| `ourgeoef` | `hypergeo_from_tables` | Complete |

## Infrastructure not compiled

- Rcpp registration and R list conversion
- Armadillo containers
- OpenMP scheduling
- R documentation-generation infrastructure
- bundled `.rda` data loading

The original files are retained under `original/` for provenance.
