# API map

| MASS/R area | Fortran API | Notes |
|---|---|---|
| `ginv`, `Null` | `ginv`, `null_space` | Pseudoinverse/eigen null space |
| `mvrnorm` | `mvrnorm` | Optional empirical moments and seed |
| `kde2d` | `kde2d` | Typed grid/density result |
| `area`, `fbeta` | `adaptive_area`, `beta_kernel` | Callback integration and beta kernel |
| `fractions`, `rational` | `rational_approximation`, `rational_values` | Integer numerator/denominator output |
| `contr.sdif` | `successive_difference_contrasts` | Exact contrast coefficients |
| `con2tr` | `con2tr` | Three-way contingency conversion |
| `lm.gls` | `lm_gls` | Explicit covariance/precision matrix |
| `lm.ridge`, `select` | `lm_ridge`, `ridge_coefficients`; inspect `gcv` | Typed ridge path |
| `rlm` | `rlm_fit` | Huber/Hampel/bisquare M fitting |
| `lqs`, `lmsreg`, `ltsreg` | `lqs_fit`, `lmsreg`, `ltsreg` | Random-subsample robust regression |
| `stdres`, `studres` | `standardized_residuals`, `studentized_residuals` | Typed regression result input |
| `boxcox`, `logtrans` | `boxcox_profile`, `logtrans_profile` | Profile log likelihood arrays |
| `dose.p` | `dose_p` | Logit/probit/cloglog |
| `huber`, `hubers` | `huber_location`, `huber_location_scale` | Typed location/scale result |
| `cov.trob` | `cov_trob` | Iteratively reweighted t covariance |
| `cov.mcd`, `cov.mve`, `cov.rob` | Same names | Vendored robust covariance foundation |
| `lda`, `qda` | `lda_fit`, `qda_fit`, prediction routines | Moment/MLE/MVE/t options where applicable |
| `fitdistr` | `fit_distribution` | Named distribution string and typed fit |
| `rnegbin` | `rnegbin` | Gamma-Poisson mixture |
| `theta.ml`, `theta.mm`, `theta.md` | `theta_ml`, `theta_mm`, `theta_md` | Explicit response and fitted means |
| `gamma.shape` | `gamma_shape_estimate` | Returns shape and standard error |
| `glm.nb` | `glm_nb_fit` | Log-link IRLS plus theta iteration |
| `negative.binomial`, `neg.bin` | `negative_binomial_variance`, `negative_binomial_deviance`, `negative_binomial_logpmf` | Numeric family components |
| `loglm`, `loglm1` | `loglinear_fit` | Explicit model matrix |
| `polr` | `polr_fit`, `polr_predict` | Five cumulative links |
| `corresp` | `correspondence_analysis` | Typed scores/correlations |
| `mca` | `mca_fit`, prediction routines | Integer factor codes and levels |
| `isoMDS` | `iso_mds` | Monotone disparities and gradient descent |
| `sammon` | `sammon` | Sammon stress minimization |
| `Shepard` | `shepard` | Target/fitted/disparity vectors |
| `ucv`, `bcv`, `width.SJ` | `ucv_bandwidth`, `bcv_bandwidth`, `width_sj` | Direct pairwise kernels |
| `bandwidth.nrd`, `nclass.freq` | `bandwidth_nrd`, `nclass_freq` | Numeric utilities |
| `addterm`, `dropterm`, `stepAIC` | `addterm_linear`, `dropterm_linear`, `step_aic_linear` | Matrix-column equivalent; no formulas |
| `negexp.SSival` | `negexp_initial` | Returns `[b0,b1,theta]` |
| graphics and S3/formula utilities | omitted | See `TRANSLATION_COVERAGE.md` |
