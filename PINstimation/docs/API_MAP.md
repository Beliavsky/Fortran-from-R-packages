# API map

## Original PIN

| R interface | Fortran counterpart | Coverage |
|---|---|---|
| `fact_pin_e()` | `pin_loglik_e()` | Direct stable factorization |
| `fact_pin_lk()` | `pin_loglik_lk()` | Direct stable factorization |
| `fact_pin_eho()` | `pin_loglik_eho()` | Direct factorization with stable fallback |
| original likelihood | `pin_loglik()` | Direct log-mixture form |
| `pin()` | `fit_pin()` | Multi-start constrained MLE |
| `pin_ea()` | `pin_ea()` | Adapted EA initializer followed by MLE |
| `pin_gwj()` | `pin_gwj()` | Adapted GWJ moment/quantile initializer followed by MLE |
| `pin_yz()` | `pin_yz()` | Adapted grid-start path followed by MLE |
| `pin_bayes()` | `fit_pin_bayes()` | Random-walk Metropolis implementation; not a bitwise port of the R/coda workflow |
| PIN posteriors | `pin_posteriors()` | Direct three-state posterior probabilities |

## Multilayer PIN

| R interface | Fortran counterpart | Coverage |
|---|---|---|
| `fact_mpin()` | `mpin_loglik()` | Direct Ersan factorization |
| `mpin_ml()` | `fit_mpin_ml()` | Constrained multi-start MLE |
| `mpin_ecm()` | `fit_mpin_ecm()` | ECM with exact posterior/weight updates and numerical conditional rate M-step |
| `get_posteriors()` | `mpin_posteriors()` | Direct no-information/good/bad layer posteriors |
| `generatedata_mpin()` | `simulate_mpin()` | Direct mixture-Poisson simulation |
| `detectlayers_e()` | `detectlayers_e()` | Adapted one-dimensional gap detector |
| `detectlayers_eg()` | `detectlayers_eg()` | Adapted one-dimensional k-means/BIC detector |
| `detectlayers_ecm()` | `detectlayers_ecm()` | BIC selection over fitted ECM models |

## Adjusted PIN

| R interface | Fortran counterpart | Coverage |
|---|---|---|
| `fact_adjpin()` | `adjpin_loglik()` | Direct stable six-component mixture likelihood |
| `adjpin(method="ML")` | `fit_adjpin_ml()` | Restricted or unrestricted constrained MLE |
| `adjpin(method="ECM")` | `fit_adjpin_ecm()` | Six-state ECM with numerical complete-data maximization |
| `generatedata_adjpin()` | `simulate_adjpin()` | Direct six-state Poisson simulation |
| AdjPIN/PSOS formulas | `adjpin_values()` | Direct formulas |
| cluster distribution | `adjpin_distribution()` | Direct six-state probabilities |

## VPIN, IVPIN, and trade data

| R interface | Fortran counterpart | Coverage |
|---|---|---|
| `classify_trades()` | `classify_trades()` | Tick, Quote, LR, and EMO algorithms; numeric timestamps |
| `aggregate_trades()` | `aggregate_classifications()` | Aggregation by caller-supplied integer group identifiers |
| VPIN bucket construction | `build_volume_buckets()` | Volume-conserving time-bar splitting and bulk classification |
| `vpin()` | `compute_vpin()` / `compute_vpin_from_buckets()` | Rolling VPIN |
| IVPIN factorization | `ivpin_loglik()` | Direct Ke-Lin likelihood convention used upstream |
| `ivpin()` | `compute_ivpin_from_buckets()` | Rolling constrained MLE with warm starts |

## Omitted R infrastructure

S4 output classes and methods, data frames and formulas, display formatting,
progress bars, parallel `future`/`furrr` execution, packaged `.RData` objects,
website/vignette generation, and calendar-aware raw-tick time-bar completion are
not native Fortran interfaces. The full upstream implementation remains in the
retained snapshot.
