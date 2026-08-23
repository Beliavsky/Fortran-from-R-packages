# API map

This file maps the important upstream `greybox` API to the modern Fortran translation.

| Upstream | Fortran | Status |
|---|---|---|
| `dlaplace/plaplace/qlaplace/rlaplace` | same lower-case names | implemented |
| `dalaplace/palaplace/qalaplace/ralaplace` | same | implemented |
| `dgnorm/pgnorm/qgnorm/rgnorm` | same | implemented |
| `ds/ps/qs/rs` | same | implemented |
| `dfnorm/pfnorm/qfnorm/rfnorm` | same | implemented |
| `dbcnorm/pbcnorm/qbcnorm/rbcnorm` | same | implemented |
| `dlogitnorm/plogitnorm/qlogitnorm/rlogitnorm` | same | implemented |
| `drectnorm/prectnorm/qrectnorm/rrectnorm` | same | implemented |
| `dtplnorm/ptplnorm/qtplnorm/rtplnorm` | same | implemented |
| `ME/MAE/MSE/MRE/MIS/MPE/MAPE/MASE/RMSSE` | lower-case equivalents | implemented |
| `rMAE/rRMSE/rAME/rMIS` | lower-case equivalents | implemented |
| `sMSE/sPIS/sCE/sMIS/GMRAE/SAME` | lower-case equivalents | implemented |
| `pinball` | `pinball` | implemented |
| `hm/ham/asymmetry/extremity/cextremity` | same | implemented |
| `measures` | `measures` | implemented numerical vector |
| `polyprod` | `polyprod` | implemented |
| `B` | `backshift` | implemented |
| `xregExpander` | `xreg_expander` | implemented numeric lag/lead expansion |
| `xregMultiplier` | `xreg_multiplier` | implemented |
| `xregTransformer` | `xreg_transformer` | implemented |
| `multipliers` | `dyn_mult_calc` | implemented recursion core |
| `outlierdummy` | `outlier_dummy` | implemented |
| `temporaldummy` | `temporal_dummy` | implemented numeric periodic dummies |
| `cramer` | `cramer_v` | implemented |
| `mcor` | `mcor` | implemented |
| `pcor` | `pcor_matrix` | implemented |
| `determination` | `determination` | implemented |
| numeric `association` | `association_numeric` | implemented |
| `alm` | `alm_fit` / `alm_model` | implemented matrix-native numerical core |
| `alm(..., distribution="dbeta")` | `alm_fit(...,'dbeta',...)` | implemented two predictor blocks |
| `alm(..., occurrence=...)` | `alm_fit_occurrence` / `alm_occurrence_model` | implemented composite occurrence-positive model |
| `alm(..., orders=c(p,d,q))` | `alm_fit_arima_errors` / `alm_dynamic_model` | implemented conditional ARIMA-error iteration |
| `alm` LASSO/RIDGE | `alm_fit(...,loss='LASSO'/'RIDGE')` | implemented |
| `alm` ROLE/QUALE | `alm_fit(...,loss='ROLE'/'QUALE')` | implemented |
| `predict.alm` | `alm_predict` | implemented |
| `AICc/BICc` | fields of `alm_model` / `alm_information` | implemented |
| `pointLik` | `alm_model%point_loglik` | implemented |
| `pAIC/pAICc/pBIC/pBICc` | `point_aic` etc. | implemented |
| `coefbootstrap` | `coef_bootstrap` | implemented |
| scale regression core | `scale_model_fit` | implemented |
| `stepwise` | `stepwise_fit` | implemented matrix-native forward IC selection |
| `lmCombine/calm` | `calm_fit`, `calm_predict` | implemented exhaustive IC averaging, including beta regression |
| `lmDynamic` | `lm_dynamic_fit` / `calm_dynamic_model` | implemented point-IC weights and LOWESS-style smoothing |
| `ro` | `rolling_origin_alm` | implemented one-step rolling-origin core |
| dynamic regression core | `recursive_lm` | implemented RLS/forgetting-factor version |
| `rmcb` | `rmcb_test` / `rmcb_result` | implemented Tukey, normal and ALM branches |
| `dsrboot` | `dsr_bootstrap` / `dsrboot_result` | implemented additive/multiplicative, parametric/nonparametric, intermittent branches |
| `aid` | `aid_fit` / `aid_result` | implemented numerical demand classification and stockout flags |
| `aidCat` | `aid_cat` / `aidcat_result` | implemented |
| `spread/tableplot/graphmaker` | -- | plotting; intentionally excluded |
| formula/S3/zoo/texreg/xtable methods | -- | R-specific; intentionally excluded |

## Deliberate interface differences

Fortran routines accept explicit numeric matrices instead of R formulas/model frames. Beta regression stores its second predictor block in `alm_model%scale_beta`; its covariance blocks are in `vcov`, `vcov_scale`, and `vcov_cross`. Occurrence and ARIMA-error models use dedicated derived types rather than overloading one large S3 object.
