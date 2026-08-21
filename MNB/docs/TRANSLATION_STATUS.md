# Translation status

Source: R package `MNB` 1.2.0 (2025-03-05).

| Upstream routine | Fortran API | Status |
|---|---|---|
| `l.MNB` | `mnb_loglik` | Complete numerical translation |
| `fit.MNB` | `fit_mnb` | Complete numerical model; native BFGS/Hessian |
| `rMNB` | `simulate_mnb` | Complete; GLG draw simplified exactly to Gamma frailty |
| `qMNB` / `f.Ymas` | `randomized_quantile_residuals`, `nb_total_pmf` | Complete |
| `re.MNB` | `residuals_mnb` | Complete source-formula translation |
| `global.MNB` | `global_influence_mnb` | Complete computational translation |
| `local.MNB`, `cases`, `cases.obs`, `cova.pertu`, `dispersion` | `local_influence_mnb` | Complete computational translation |
| `envelope.MNB` | `envelope_mnb` | Complete computational translation |

## Omitted infrastructure

- plotting and normal-probability/index plots;
- R formula/model-frame/model-matrix parsing;
- data-frame manipulation and named-list presentation;
- bundled `.rda` datasets as Fortran constants;
- R-specific RNG stream identity and exact `optim` line-search micro-behavior.

The statistical formulas and diagnostic calculations themselves are translated.
