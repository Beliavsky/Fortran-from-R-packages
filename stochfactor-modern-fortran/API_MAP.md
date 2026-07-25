# Computational API map

This file maps important original R-facing routines to the combined Fortran API.
"Analogue" means the numerical purpose is implemented but the algorithm or full
option surface is not identical to R/C++.

## stochvol

| Original routine or family | Fortran procedure/type | Status |
|---|---|---|
| `svsim` | `simulate_sv`, `sv_sim_result` | Implemented |
| Gaussian SV | `svsample` with default options | Numerical analogue |
| `svtsample` | `svsample_t` | Numerical analogue |
| `svlsample` | `svsample_leverage` | Numerical analogue |
| `svtlsample` | `svsample_t_leverage` | Numerical analogue |
| `svsample`, `svsample2` | `svsample` | Numerical analogue |
| `svsample_fast_cpp` | Omori constants, indicator draw, `draw_latent_mixture` through `svsample` | Partial numerical analogue; no original ASIS/adaptation engine |
| `svsample_general_cpp` | general latent RW and parameter updates through `svsample` | Partial numerical analogue |
| `predict.svdraws`, `predy`, `predlatent`, `predvola` | `predict_sv`, `sv_prediction` | Implemented numerical outputs |
| `residuals.svdraws` | `sv_residuals` | Implemented |
| rolling sample wrappers | `rolling_sv_forecast` | Simplified expanding rolling analogue |
| regression-mean updates | `bayesian_regression_update`, beta updates in `svsample` | Implemented |
| `svtau`, `svbeta`, `latent`, `latent0`, `para` | fields of `sv_draws` | Implemented as plain arrays |
| `logret` | `log_returns` | Implemented |
| `specify_priors`, prior classes | `sv_prior` | Partial compact replacement |
| plotting, S3, `coda`, summaries | none | Excluded infrastructure |

## factorstochvol

| Original routine or family | Fortran procedure/type | Status |
|---|---|---|
| `fsvsim` | `simulate_fsv`, `fsv_sim_result` | Implemented |
| `fsvsample` | `fit_fsv`, `fsv_options`, `fsv_draws` | Broad numerical analogue; no interweaving parity |
| default dense/sparse simulation loadings | `make_dense_loadings`, `make_sparse_loadings` | Implemented |
| default SV parameters | `default_fsv_parameters` | Implemented |
| `covmat` | `fsv_covariance_path` | Implemented |
| `cormat` | `fsv_correlation_path` | Implemented |
| `covelement` | `covelement` | Implemented |
| `corelement` | `corelement` | Implemented |
| `expweightcov` | `expweightcov` | Implemented |
| `ledermann` | `ledermann` | Implemented |
| factor initialization | `static_factor_initialize` | PCA-based analogue |
| `preorder` | `preorder` | Numerical analogue |
| `findrestrict` | `findrestrict` | Numerical analogue |
| loading restrictions | optional `restriction` argument to `fit_fsv` | Implemented |
| Normal-Gamma prior | `fsv_options%normal_gamma`, stored local/global scales | Implemented with GIG slice sampler |
| `predh`, `predcov`, `predcor` | `predict_fsv`, `fsv_prediction` | Implemented |
| `predprecWB` | `woodbury_precision`, prediction precision fields | Implemented |
| `predloglik` | `predloglik_fsv` | Implemented |
| `predloglikWB` | `predloglik_fsv_woodbury` | Implemented |
| predictive draw aggregation | `aggregate_loglik_draws` | Implemented |
| `predcond` | `predict_conditional_fsv` | Implemented numerical outputs |
| internal `dmvnorm`/vectorized density | `dmvnorm_columns` | Implemented |
| `signident` | `sign_identify` | Implemented |
| `orderident` | `order_identify` | Implemented |
| `evdiag` | `eigen_loading_diagnostics` | Implemented |
| `runningcovmat` | `running_covariance` | Implemented |
| `runningcormat` | `running_correlation` | Implemented |
| plots, R methods, dates, parallel sampler plumbing | none | Excluded infrastructure |

## Supporting modules

- `sv_rng`: deterministic seeding, Uniform, Normal, Gamma, inverse-Gamma,
  standardized Student-t, and GIG log-scale slice generation.
- `sv_linalg`: Cholesky decomposition, SPD solve/inverse/log determinant,
  multivariate Normal density/draw, symmetric eigendecomposition, covariance and
  correlation matrices.
- `sv_stats`: basic summaries, Normal and standardized Student-t log densities,
  stable log-sum-exp, quantiles, and log returns.
