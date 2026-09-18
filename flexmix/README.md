# flexmix — modern Fortran translation

This directory translates the computational core of the R package **flexmix 2.3-21**
to modern free-form Fortran with the Fortran Package Manager (FPM).

The upstream package implements a general finite-mixture framework. This translation
focuses on its numerical EM machinery and the most reusable regression, clustering,
concomitant-prior, inference, and model-selection operations. It intentionally accepts
numeric arrays rather than reproducing R formulas, S4 classes, model frames, plotting,
or interactive interfaces.

## Implemented numerical scope

The public `flexmix` module provides:

- Gaussian, Poisson, binomial, and Gamma inverse-link mixtures of regressions, with numeric offsets.
- Multinomial-logit response mixtures with a class-one reference parameterization.
- Zero-inflated Poisson/binomial and robust-background Gaussian/Poisson GLM mixtures.
- Full/diagonal multivariate-normal clustering.
- Independent multivariate Bernoulli and Poisson mixtures.
- Mixed Bernoulli/independent-Gaussian clustering.
- Lognormal, exponential, inverse-Gaussian, Gamma, and Weibull univariate mixtures.
- Weighted or hard EM classification, grouped observations, `minprior` component
  deletion, and optional multinomial-logit concomitant priors.
- Fitted priors, posteriors, clusters, parameters, component predictions, relabeling,
  component removal, Gaussian refitting, packed parameter replacement, analytic regression
  scores, and observed-information covariance/standard errors.
- Numeric fixed/shared-coefficient Gaussian, Poisson, binomial, and Gamma mixture drivers,
  including shared Gaussian variance groups.
- Self-contained elastic-net/lasso Gaussian, Poisson, and binomial mixtures with adaptive
  penalties, offsets, selection masks, fixed-fold cross-validation, and lambda selection.
- Grouped Gaussian linear mixed-model mixtures with arbitrary numeric random-effect design
  matrices, plus a partial `FLXMRlmer` analogue for the same marginal model.
- Left-censored Gaussian mixed-model mixtures with exact one-dimensional truncation and
  deterministic multivariate truncated-normal integration for several censored rows per group.
- Penalized smooth Gaussian, Poisson, and binomial mixtures from caller-supplied basis/design
  and smoothing-penalty matrices, with selected smoothing parameters and effective df.
- `logLik`, AIC, BIC, ICL, EIC, discrete/model KL divergence, deterministic repeated
  Gaussian model selection, step-result deduplication, regression bootstrap, and adjacent-
  component bootstrap likelihood-ratio testing.
- `ExLinear`, `ExNPreg`, and `ExNclus`-style data generation plus simulation from supported
  regression fits and fixed numerical mixture distributions.

No external numerical library is needed. The package does not link BLAS/LAPACK and does
not vendor `rfortran-*`, `nnet`, `mclust`, or any other translated dependency.

## Build

With FPM and gfortran available:

```text
fpm build
fpm test
fpm run --example basic_flexmix
```

The maintained Fortran sources use a single `dp = real64` kind defined in
`flexmix_kinds`; all translated/new source is free form.

## Example

```fortran
use flexmix, only : dp, flexmix_result, flexmix_gaussian

real(dp) :: x(20,2), y(20)
type(flexmix_result) :: fit

! Fill x and y, including an intercept column in x.
call flexmix_gaussian(x, y, 2, fit)
print *, fit%prior
print *, fit%beta
```

See `example/basic_flexmix.f90` for a complete deterministic program.

## Translation coverage

Package status: **substantial**.

**66 of 68 (97.1%)** exported computational R functions/methods have a meaningful
Fortran mapping. This percentage measures mapped exported computational functions,
**not** complete R-interface compatibility. The remaining gaps are R formula/model-matrix and callback collection infrastructure.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `ExNPreg` | `flexmix_api` | `ex_npreg` | substantial |
| `ExNclus` | `flexmix_api` | `ex_nclus` | substantial |
| `FLXdist` | `flexmix_api` | `flexmix_make_distribution` | partial |
| `FLXgradlogLikfun` | `flexmix_refit` | `flexmix_gradient_regression` | substantial |
| `FLXglmFix` | `flexmix_fixed` | `flexmix_glmfix_gaussian`, `flexmix_glmfix_poisson`, `flexmix_glmfix_binomial`, `flexmix_glmfix_gamma` | substantial |
| `FLXMRglmfix` | `flexmix_fixed` | `flexmix_glmfix_gaussian`, `flexmix_glmfix_poisson`, `flexmix_glmfix_binomial`, `flexmix_glmfix_gamma` | substantial |
| `FLXgetDesign` | `flexmix_refit` | `flexmix_get_design` | substantial |
| `FLXreplaceParameters` | `flexmix_refit` | `flexmix_replace_parameters` | substantial |
| `existGradient` | `flexmix_refit` | `flexmix_exist_gradient` | substantial |
| `refit_optim` | `flexmix_refit` | `flexmix_refit_optim_regression` | partial |
| `boot` | `flexmix_bootstrap` | `flexmix_boot_regression` | partial |
| `LR_test` | `flexmix_bootstrap` | `flexmix_lr_test_regression` | partial |
| `unique` | `flexmix_api` | `step_flexmix_unique` | substantial |
| `ExLinear` | `flexmix_api` | `ex_linear` | partial |
| `FLXlogLikfun` | `flexmix_api`, `flexmix_refit` | `flexmix_loglik`, `flexmix_loglik_regression_at` | substantial |
| `FLXMCfactanal` | `flexmix_api` | `flexmix_factanal` | substantial |
| `FLXMCdist1` | `flexmix_api` | `flexmix_dist1` | substantial |
| `FLXMCmvbinary` | `flexmix_api` | `flexmix_mvbinary` | substantial |
| `FLXMCmvcombi` | `flexmix_api` | `flexmix_mvcombi` | substantial |
| `FLXMCmvnorm` | `flexmix_api` | `flexmix_mvnorm` | substantial |
| `FLXMCmvpois` | `flexmix_api` | `flexmix_mvpois` | substantial |
| `FLXMCnorm1` | `flexmix_api` | `flexmix_norm1` | substantial |
| `FLXPmultinom` | `flexmix_em_utils` | `fit_multinomial_prior`, `multinomial_prior` | substantial |
| `FLXPconstant` | `flexmix_em_utils` | `fit_constant_prior` | substantial |
| `FLXfit` | `flexmix_core` | `flexmix_fit_regression`, `flexmix_fit_binomial`, `flexmix_fit_multinomial`, `flexmix_fit_conditional_logit`, `flexmix_fit_multivariate`, `flexmix_fit_univariate` | substantial |
| `FLXglm` | `flexmix_api` | `flexmix_gaussian`, `flexmix_poisson`, `flexmix_binomial`, `flexmix_gamma` | substantial |
| `FLXbclust` | `flexmix_api` | `flexmix_mvbinary` | substantial |
| `FLXmclust` | `flexmix_api` | `flexmix_mvnorm` | substantial |
| `FLXconstant` | `flexmix_em_utils` | `fit_constant_prior` | substantial |
| `FLXmultinom` | `flexmix_em_utils` | `fit_multinomial_prior`, `multinomial_prior` | substantial |
| `FLXMRcondlogit` | `flexmix_api` | `flexmix_conditional_logit` | substantial |
| `FLXMRglm` | `flexmix_api` | `flexmix_gaussian`, `flexmix_poisson`, `flexmix_binomial`, `flexmix_gamma` | substantial |
| `FLXMRglmnet` | `flexmix_api` | `flexmix_glmnet_gaussian`, `flexmix_glmnet_poisson`, `flexmix_glmnet_binomial` | substantial |
| `FLXMRlmm` | `flexmix_mixed` | `flexmix_lmm` | substantial |
| `FLXMRlmmc` | `flexmix_mixed` | `flexmix_lmc`, `flexmix_lmmc` | substantial |
| `FLXMRlmer` | `flexmix_mixed` | `flexmix_lmer` | partial |
| `FLXMRmgcv` | `flexmix_api` | `flexmix_mgcv_gaussian`, `flexmix_mgcv_poisson`, `flexmix_mgcv_binomial` | partial |
| `FLXMRmultinom` | `flexmix_api` | `flexmix_multinomial` | substantial |
| `FLXMRrobglm` | `flexmix_api` | `flexmix_robust_gaussian`, `flexmix_robust_poisson` | substantial |
| `FLXMRziglm` | `flexmix_api` | `flexmix_ziglm_poisson`, `flexmix_ziglm_binomial` | substantial |
| `FLXgetParameters` | `flexmix_api` | `flexmix_get_parameters` | substantial |
| `flexmix` | `flexmix_api` | `flexmix_gaussian`, `flexmix_poisson`, `flexmix_binomial`, `flexmix_gamma`, `flexmix_multinomial`, `flexmix_conditional_logit`, `flexmix_ziglm_poisson`, `flexmix_ziglm_binomial`, `flexmix_robust_gaussian`, `flexmix_robust_poisson`, `flexmix_mvnorm`, `flexmix_factanal`, `flexmix_mvbinary`, `flexmix_mvpois`, `flexmix_mvcombi`, `flexmix_dist1` | substantial |
| `group` | `flexmix_em_utils` | `group_first_mask`, `group_log_densities` | partial |
| `initFlexmix` | `flexmix_api` | `step_flexmix_gaussian` | partial |
| `relabel` | `flexmix_api` | `flexmix_relabel` | complete |
| `stepFlexmix` | `flexmix_api` | `step_flexmix_gaussian` | partial |
| `EIC` | `flexmix_api` | `flexmix_eic` | complete |
| `FLXcheckComponent` | `flexmix_api` | `flexmix_check_result` | partial |
| `FLXdeterminePostunscaled` | `flexmix_components` | `gaussian_regression_log_density`, `poisson_regression_log_density`, `binomial_regression_log_density`, `gamma_regression_log_density`, `multinomial_regression_log_density`, `conditional_logit_log_density`, `mvnormal_log_density`, `mvbinary_log_density`, `mvpois_log_density`, `mvcombi_log_density`, `lognormal_log_density`, `exponential_log_density`, `inverse_gaussian_log_density`, `gamma_log_density`, `weibull_log_density` | substantial |
| `FLXgetK` | `flexmix_api` | `flexmix_get_k` | complete |
| `FLXgetObs` | `flexmix_api` | `flexmix_get_obs` | complete |
| `FLXmstep` | `flexmix_components` | `fit_gaussian_regression_component`, `fit_poisson_regression_component`, `fit_binomial_regression_component`, `fit_gamma_regression_component`, `fit_multinomial_regression_component`, `fit_conditional_logit_component`, `fit_mvnormal_component`, `fit_factor_analysis_component`, `fit_mvbinary_component`, `fit_mvpois_component`, `fit_mvcombi_component`, `fit_lognormal_component`, `fit_exponential_component`, `fit_inverse_gaussian_component`, `fit_gamma_component`, `fit_weibull_component` | substantial |
| `FLXremoveComponent` | `flexmix_api` | `flexmix_remove_component` | partial |
| `KLdiv` | `flexmix_api` | `kl_divergence_matrix`, `kl_divergence_regression`, `kl_divergence_mvnorm` | substantial |
| `ICL` | `flexmix_api` | `flexmix_icl` | complete |
| `clusters` | `flexmix_api` | `flexmix_get_clusters` | substantial |
| `fitted` | `flexmix_api` | `flexmix_predict_regression`, `flexmix_predict_multinomial` | partial |
| `getModel` | `flexmix_api` | `get_model` | substantial |
| `logLik` | `flexmix_api` | `flexmix_loglik` | substantial |
| `parameters` | `flexmix_api` | `flexmix_get_parameters` | substantial |
| `posterior` | `flexmix_api` | `flexmix_get_posterior` | partial |
| `predict` | `flexmix_api` | `flexmix_predict_regression`, `flexmix_predict_multinomial` | partial |
| `prior` | `flexmix_api` | `flexmix_get_prior` | substantial |
| `rFLXM` | `flexmix_api` | `flexmix_simulate_regression` | partial |
| `refit` | `flexmix_api` | `flexmix_refit_gaussian` | partial |
| `rflexmix` | `flexmix_api` | `flexmix_simulate_regression` | partial |

The 2 functions still outside the translated computational surface are:
`FLXgetModelmatrix` and `Lapply`.

The machine-readable mapping in `fpm.toml` is authoritative and uses the same counts.

## Numerical and interface differences

- Initial partitions are deterministic unless the caller supplies `initial_cluster` or
  `initial_posterior`; upstream R may use random starts depending on the workflow.
- Numeric design matrices are passed directly. The caller must create intercept,
  factor/dummy, and transformed columns that R formulas would otherwise create; scalar
  GLM offsets are accepted explicitly by the corresponding Fortran entry points.
- Gaussian weighted fitting follows the upstream `lm.wfit` residual-scale convention
  used by `FLXMRglm`. Multivariate normal covariance follows the upstream
  `cov.wt(..., method="unbiased")` semantics.
- Both multinomial concomitant priors and multinomial response mixtures use internal
  Newton/Fisher-scoring solvers rather than `nnet::multinom`.
- `step_flexmix_gaussian` currently covers Gaussian regression mixtures only;
  `step_flexmix_unique` operates on that numeric step-result representation.
- `flexmix_boot_regression`/`flexmix_lr_test_regression` cover ordinary regression
  mixtures without grouped observations or concomitant priors. Parametric Gamma
  bootstrap generation uses a Marsaglia-Tsang sampler.
- `flexmix_simulate_regression` and the example/bootstrap generators do not attempt to
  reproduce R's RNG stream. Poisson simulation uses exact Knuth sampling at moderate
  means and a normal approximation above 30.
- The fixed/shared-coefficient GLM API uses explicit component-specific numeric design
  arrays and an incidence matrix instead of nested R formulas; automatic component
  deletion is intentionally disabled for that fixed design.
- `FLXMRglmnet` uses an internal weighted proximal-gradient elastic-net solver. It does
  not reproduce `glmnet` standardization/path generation or R's random CV-fold defaults.
- `FLXMRlmm` translates the grouped Gaussian random-effect EM updates on explicit fixed
  and random design matrices. `FLXMRlmer` exposes the same marginal model numerically,
  but does not reproduce `lme4` formula syntax or optimizer internals.
- `FLXMRlmmc` now covers both upstream numerical branches: `flexmix_lmc` handles censored
  Gaussian regression without random effects, while `flexmix_lmmc` handles grouped random
  effects. A single censored row in a mixed-model group uses exact conditional truncated-normal
  moments; several censored rows use a deterministic Genz-style Halton QMC integral. The
  `mvtnorm` integration-control surface and R formula/model-frame construction are not reproduced.
- `FLXMRmgcv` accepts a preconstructed numeric basis/design matrix and one symmetric
  smoothing-penalty matrix. It does not construct `mgcv::s()` bases or reproduce the
  full multi-penalty `mgcv` smoothing-parameter optimizer.

## Provenance and license

See `NOTICE.md`, `LICENSE.md`, and `COPYING`. Original R sources needed to audit the
function mapping are preserved under `upstream/R/` with their original notices.
