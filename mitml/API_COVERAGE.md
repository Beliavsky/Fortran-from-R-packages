# API coverage

This document maps the computational surface of upstream **mitml 0.4-5** to
the maintained Fortran API. Plotting, file I/O, S3 dispatch, formula parsing,
model-object extraction/refitting, and dynamic R expression evaluation are
intentionally excluded.

| Upstream computation | Fortran API | Coverage |
|---|---|---|
| `.pool.estimates` / numeric `testEstimates` | `pool_estimates`, `pool_confint` | Implemented, including total covariance, RIV, FMI, t tests, Barnard-Rubin finite complete-data df |
| `.D1` | `d1_test` | Implemented, including Reiter finite complete-data df correction |
| `.D2` | `d2_test` | Implemented |
| `testModels(..., method="D1")` numerical combination | `d1_test` | Implemented once estimates/covariances have been extracted from fits |
| `testModels(..., method="D2")` | `d2_test` | Implemented for supplied per-imputation Wald or likelihood-deviance values |
| `testModels(..., method="D3")` | `d3_test` | Implemented for supplied per-imputation and pooled-parameter log likelihoods |
| `testModels(..., method="D4")` | `d4_test` | Implemented, including default/positive RIV and robust Chan-Meng RIV paths |
| `testConstraints` linear restrictions | `test_linear_constraints` | Implemented |
| `testConstraints` arbitrary parsed R expressions | `test_transformed_constraints` | Numerical pooling implemented after caller supplies transformed estimates/Jacobian covariance; R expression parser and `numericDeriv` orchestration intentionally omitted |
| `clusterMeans` | `cluster_means`, `cluster_means_matrix` | Implemented, including NaNs, nested `group`, and `adj=TRUE` leave-one-out means |
| `.getRsquared` numerical formulas used by `multilevelR2` | `multilevel_r2`, `intraclass_correlation` | Implemented for supplied model components |
| Null-model refitting and `lme4`/`nlme` extraction in `multilevelR2` | — | R/model-interface specific; omitted |
| `.GelmanRubin` | `gelman_rubin` | Implemented with upstream small-sample variance and df calculation |
| `.SDprop` | `sd_proportion` | Implemented using the same spectral-density-at-zero definition and automatic Yule-Walker AIC order selection; tiny order-selection details may differ from `stats::ar` |
| `.reducedACF` / `.smoothedACF` numerical kernel | `reduced_acf` | Implemented with normal-kernel lag smoothing |
| `.movingAverage` | `moving_average` | Implemented, including upstream edge-fill behavior |
| `.logLik_lm` | `gaussian_lm_loglik` | Implemented exactly, including the Gaussian normalizing constant |
| `.logLik_lmm` | `gaussian_lmm_loglik` | Implemented for supplied numeric design matrices, cluster labels, fixed effects, random-effect covariance, and residual variance; matches upstream helper, including omission of the constant `-n*log(2*pi)/2` |
| `panImpute` | separate top-level `pan` translation | Upstream function is R formula/data orchestration around `pan`; not linked here |
| `jomoImpute` | separate top-level `jomo` translation | Upstream function is R formula/data orchestration around `jomo`; not linked here |
| `mitmlComplete`, `as.mitml.list`, converters, `with`/`within`, row/column binding | — | R list/data-frame management; omitted |
| `read.mitml`, `write.mitml*` | — | R/file-format interface; omitted |
| plot/print/summary methods | — | Plotting, presentation, and S3 interface code omitted; reusable convergence numerics retained |
| model extraction helpers (`.getCoef`, `.getVcov`, `.getLL`, `.updateML`, etc.) | — | R class dispatch/model-object plumbing; callers supply numerical arrays directly |

## Numerical design differences

Fortran has no safe equivalent of evaluating arbitrary user-supplied R
expressions by `parse()` and `numericDeriv()`. General nonlinear constraints
are therefore split into two explicit stages: callers evaluate the constraints
and their Jacobians in application code, then pass the transformed estimates
and covariance matrices to `test_transformed_constraints`.

`sd_proportion` fits Yule-Walker AR candidates and chooses the order by an AIC
criterion. This preserves mitml's statistical definition of the diagnostic,
but is not intended to reproduce every internal detail of the R `stats::ar`
implementation bit-for-bit.
