# Translation coverage

Upstream version: `glmnet 5.0`.

| Upstream export | Fortran representation | Status |
|---|---|---|
| `glmnet` | `fit_glmnet`, specialized family fits | Implemented |
| `cv.glmnet` | `cv_glmnet`, `cv_multinomial`, `cv_cox` | Implemented |
| `predict.glmnet` | `predict_glmnet`, `predict_glmnet_at` | Implemented |
| `coef.glmnet` | `coef_glmnet` | Implemented |
| `relax.glmnet` | four typed relax procedures | Implemented |
| `predict.relaxed` | `predict_glmnet` on relaxed result | Implemented |
| `coef.relaxed` | `coef_glmnet` on relaxed result | Implemented |
| `assess.glmnet` | three typed assessment procedures | Implemented |
| `confusion.glmnet` | `confusion_glmnet` | Implemented |
| `roc.glmnet` | `roc_glmnet` | Implemented |
| `Cindex` | `concordance_index` | Implemented |
| `coxgrad` | `cox_gradient` | Implemented |
| `coxnet.deviance` | `coxnet_deviance` | Implemented |
| `glmnet.control` | `glmnet_control_type`, update helpers | Implemented |
| `glmnet.measures` | `glmnet_measures` | Implemented |
| `bigGlm` | `big_glm` | Implemented |
| `buildPredmat` | `build_predmat`, CV result predictions | Implemented |
| `makeX` | `make_x` | Numeric-matrix equivalent |
| `prepareX` | `prepare_x` | Numeric-matrix equivalent |
| `na.replace` | `na_replace` | Implemented |
| `na_sparse_fix` | `na_sparse_fix` | Implemented |
| `rmult` | `rmult` | Implemented |
| `stratifySurv` | `stratify_surv` | Implemented |
| `print.cv.glmnet` | ordinary Fortran output by caller | Display-only omission |

## Family coverage

| Family/feature | Coverage |
|---|---|
| Gaussian elastic net | Yes |
| Binomial logistic elastic net | Yes |
| Poisson elastic net | Yes |
| Multinomial, ungrouped | Yes |
| Multinomial, grouped | Yes |
| Multiresponse Gaussian grouped penalty | Yes |
| Cox right-censored | Yes |
| Cox `(start, stop]` | Yes |
| Cox strata | Yes |
| Cox Efron/Breslow ties | Yes |
| Observation weights | Yes |
| Offsets | Yes |
| Penalty factors | Yes |
| Coefficient bounds | Yes |
| Excluded predictors | Yes |
| Automatic/user lambda path | Yes |
| Dense input | Yes |
| CSC input | Dense fallback |
| Arbitrary GLM family | Typed IRLS callback |
| Relaxed path | Yes |
| Cross-validation | Sequential |

The port does not claim the upstream C++ engine's strong-rule screening,
method-specific sparse complexity, or bit-for-bit iteration path.
