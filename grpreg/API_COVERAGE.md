# API coverage

Coverage follows the repository convention: distinct exported computational R functions and explicitly exported computational S3 methods are counted. Plotting, printing, presentation-only methods, datasets, and dispatch-only generics are excluded.

**Status: substantial — 19 of 19 (100%) mapped computational functions/methods.**

| R function | Fortran module | Public procedure(s) | Status |
|---|---|---|---|
| `grpreg` | `grpreg_api` | `grpreg_fit` | substantial |
| `grpsurv` | `grpreg_api` | `grpsurv_fit` | substantial |
| `gBridge` | `grpreg_api` | `gbridge_fit` | substantial |
| `cv.grpreg` | `grpreg_cv` | `cv_grpreg` | substantial |
| `cv.grpsurv` | `grpreg_cv` | `cv_grpsurv` | substantial |
| `expand_spline` | `grpreg_spline` | `expand_spline`, `predict_spline` | substantial |
| `gen_nonlinear_data` | `grpreg_data` | `gen_nonlinear_data` | substantial |
| `mfdr` | `grpreg_mfdr` | `mfdr_grpreg` | substantial |
| `AUC.cv.grpsurv` | `grpreg_cv` | `auc_cv_grpsurv` | partial |
| `coef.grpreg` | `grpreg_api` | `coef_grpreg` | complete |
| `coef.cv.grpreg` | `grpreg_cv` | `coef_cv_grpreg` | complete |
| `logLik.grpreg` | `grpreg_api` | `loglik_grpreg` | substantial |
| `logLik.grpsurv` | `grpreg_api` | `loglik_grpreg` | complete |
| `predict.grpreg` | `grpreg_api` | `predict_grpreg`, `coef_grpreg`, `count_nonzero`, `count_nonzero_groups`, `group_norms` | substantial |
| `predict.cv.grpreg` | `grpreg_cv` | `predict_cv_grpreg` | substantial |
| `predict.grpsurv` | `grpreg_api` | `predict_grpsurv_link`, `predict_grpsurv_survival`, `predict_grpsurv_hazard`, `predict_grpsurv_median` | substantial |
| `residuals.grpreg` | `grpreg_api` | `residuals_grpreg` | substantial |
| `select.grpreg` | `grpreg_api` | `select_grpreg` | substantial |
| `summary.cv.grpreg` | `grpreg_cv` | `summarize_cv_grpreg` | substantial |

The 100% fraction means that every function in this **coverage basis** has a meaningful mapping. It does not mean complete compatibility with all R call signatures, S3 object behavior, attributes, multitask matrix responses, screening accelerators, or RNG streams.
