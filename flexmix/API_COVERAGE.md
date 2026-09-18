# API coverage

This document summarizes numerical scope beyond the per-function mapping in `fpm.toml`.
The audited exported-computational coverage is **66 of 68 (97.1%)**. That fraction records
meaningful R-to-Fortran function mappings; it does not imply complete R/S4 compatibility.

## Core EM

`flexmix_core` implements translated E/M loops for scalar regression, binomial and
multinomial regression, robust/structural-zero GLM variants, multivariate clustering,
and positive univariate distribution mixtures. It supports case weights, grouped
observations, hard or weighted classification, constant or multinomial concomitant
priors, automatic deletion of components below `minprior`, stable log-sum-exp
normalization, and the upstream relative log-likelihood convergence test.

`flexmix_fixed` adds a numeric `FLXMRglmfix`/`FLXglmFix` analogue. Callers provide one
component-specific design matrix per component plus a logical incidence matrix identifying
shared and component-specific coefficients. Gaussian variance groups can also be shared.
The driver covers Gaussian, Poisson, binomial, and Gamma inverse-link families. It does not
attempt R nested-formula construction and deliberately does not delete components because
that would invalidate the fixed component/design incidence supplied by the caller.

`flexmix_penalized` adds weighted elastic-net/lasso component fitting for Gaussian,
Poisson, and binomial regression mixtures, including adaptive penalty factors, explicit
offsets, selection masks, fixed-fold cross-validation, and per-component lambda selection.
The implementation is self-contained and does not depend on `glmnet`.

`flexmix_mixed` implements grouped Gaussian random-effect mixtures on explicit fixed and
random design matrices. `flexmix_lmm` follows the upstream EM covariance updates;
`flexmix_lmer` exposes the same marginal Gaussian model as a documented partial mapping
without reproducing `lme4` formula/optimizer machinery. For `FLXMRlmmc`, `flexmix_lmc`
implements the no-random-effects censored Gaussian branch with exact univariate truncated-
normal moments, while `flexmix_lmmc` handles random effects. In the mixed branch,
one-dimensional truncation is analytic and multiple censored rows use a deterministic
Genz/Halton QMC truncated-normal probability-and-moment kernel (up to 32 censored dimensions).

`flexmix_smooth` implements the post-formula numerical core of `FLXMRmgcv` for Gaussian,
Poisson, and binomial families. Callers supply a basis/design matrix and symmetric
positive-semidefinite smoothing penalty. The implementation returns selected smoothing
parameters and effective degrees of freedom.

## Refit and inference machinery

`flexmix_refit` provides an unconstrained numerical parameter vector, component-to-unique-
parameter design matrices, parameter replacement, arbitrary-parameter regression-mixture
log likelihoods, analytic observed-data scores for ordinary Gaussian/Poisson/binomial
mixtures, and optimizer-style observed-information covariance/standard errors. Fixed/shared
models participate in parameter design, packing, and replacement, but their analytic score
path is intentionally reported unavailable.

`flexmix_bootstrap` provides parametric and row-resampling bootstrap fits for ordinary
Gaussian, Poisson, binomial, and Gamma regression mixtures plus an adjacent-component
parametric-bootstrap likelihood-ratio test. Grouped observations, concomitant priors, and
non-regression component families remain outside this bootstrap layer.

## Supported component families

| Family | M-step | log density | top-level fit |
| --- | --- | --- | --- |
| Gaussian regression | yes | yes | yes |
| Poisson regression | yes | yes | yes |
| Binomial regression | yes | yes | yes |
| Gamma regression (inverse link) | yes | yes | yes |
| Penalized Gaussian/Poisson/binomial regression | yes | yes | yes |
| Penalized smooth Gaussian/Poisson/binomial regression | yes | yes | yes |
| Grouped Gaussian linear mixed model | yes | yes (marginal Gaussian) | yes |
| Left-censored Gaussian regression / mixed model | yes | yes (censored Gaussian) | yes |
| Fixed/shared-coefficient GLM regression | yes | yes | yes |
| Multinomial logistic regression | yes | yes | yes |
| Zero-inflated Poisson/binomial | yes | yes | yes |
| Robust-background Gaussian/Poisson | yes | yes | yes |
| Multivariate normal | yes | yes | yes |
| Factor analyzer | yes | yes (multivariate normal covariance) | yes |
| Conditional logistic regression | yes | yes | yes |
| Multivariate binary | yes | yes | yes |
| Multivariate Poisson | yes | yes | yes |
| Mixed binary/Gaussian | yes | yes | yes |
| Lognormal | yes | yes | yes |
| Exponential | yes | yes | yes |
| Inverse Gaussian | yes | yes | yes |
| Gamma | yes | yes | yes |
| Weibull | yes | yes | yes |


## Factor-analysis and conditional-choice models

`flexmix_factanal` fits mixtures of factor analyzers from numeric multivariate data. Its component M-step reproduces the `cov.wt` covariance preprocessing and the concentrated maximum-likelihood `factanal.fit.mle` uniqueness objective with bounds `[0.005,1]`; the returned loadings are unrotated and the covariance used by the likelihood is rotation invariant.

`flexmix_conditional_logit` fits one-choice-per-stratum conditional-logit mixtures from an explicit no-intercept design matrix and integer stratum labels. The M-step uses a stable stratum softmax likelihood and Newton updates, while the E-step reproduces the upstream per-row divided stratum log likelihood before grouped mixture normalization. R formula preprocessing and the more general `survival::coxph` exact-tie machinery are not reproduced.

## Major upstream features not yet translated

Two audited exports remain unmapped:

- `FLXgetModelmatrix` — R formula/model-frame construction; and
- `Lapply` — an R callback/list convenience method over refitted component objects.

Both remaining exports are R interface infrastructure. They are intentionally not claimed
merely because numeric Fortran arrays and loops can perform related lower-level operations.
