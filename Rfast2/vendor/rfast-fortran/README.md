# Rfast-fortran

Modern Fortran 2018 / FPM computational-core translation of **Rfast 2.1.5.2**.

Rfast is library-scale: the supplied NAMESPACE exports 444 unique names spanning array kernels, distributions, fitted models, hypothesis tests, model selection, R object utilities, and specialized workflows. Version 0.3.0 expands the reusable Fortran core substantially while continuing to document untranslated exports explicitly instead of claiming false full parity.

`API_MAP.md` accounts for every upstream export. In v0.3.0, **257** are ported/consolidated, **1** is partially ported, **13** map directly to Fortran intrinsics, **28** are R-only omissions, and **145** computational exports remain pending.

## v0.3.0 additions

- completed every response-family branch of upstream `omp`: logistic, Poisson, quasi-Poisson, quasi-binomial, normal-log, Gamma, Weibull, multivariate Gaussian, and multinomial
- random-intercept MLE/regression for balanced or unbalanced groups, random effects, batched single-predictor tests, and the balanced `rint.regbx` shortcut
- normal-log and Weibull regression plus their many-predictor wrappers
- Tobit MLE and ordinal logit/probit MLE
- batched Gamma, inverse-Gaussian, quasi-Poisson, proportion, and multinomial regression tests, plus numeric-matrix `univglms`
- native PC-skeleton implementation for Pearson and Spearman correlation with arbitrary analytic conditioning order
- four new v0.3 regression suites with independent reference checks for Tobit, normal-log, and Weibull likelihoods

The PC-skeleton categorical and permutation-test modes remain the single partial export in the coverage map.

## v0.2.0 additions

- all ten upstream score-test families: beta, exponential, gamma, geometric, generic GLM, inverse-Gaussian, multinomial, negative-binomial, Weibull, and zero-truncated Poisson
- multivariate normal, lognormal, Student-t, Dirichlet, and inverse-Dirichlet MLEs
- one- and two-sample empirical likelihood and exponential empirical likelihood
- multivariate exponential empirical likelihood
- repeated-measures ANOVA, batched repeated-measures ANOVA, and repeated-measures line tests
- method-of-moments variance components
- balanced and unbalanced random-intercept variance-component MLEs, including random effects
- Gamma, inverse-Gaussian, quasi-Poisson, proportion, multinomial, and spatial-median regressions
- OMP/BIC linear-model selection, partial-correlation BIC forward selection, and logistic/Poisson OMP selection

The unbalanced variance-component MLE preserves Rfast's `gold_rat3` bounded variance-ratio search rather than silently substituting a different mixed-model optimizer.

## Existing numerical core

The v0.1 functionality remains available: row/column descriptive statistics and ordering, matrix algebra, covariance/correlation, distances and energy statistics, special functions, multivariate densities/RNGs, broad scalar and column-wise MLEs, linear/logistic/Poisson/ridge regression, AR(1), classical tests and ANOVA, k-NN, graph/permutation algorithms, spatial median, directional methods, and Gaussian/Poisson/geometric/multinomial/Gamma naive Bayes.

The package vendors the previously translated `zigg-fortran` package and uses its 32-bit-compatible Ziggurat/KISS generator for normal, exponential, and uniform streams.

## Build

```text
fpm build
fpm test
fpm run --example basic_stats
fpm run --example v02_features
fpm run --example v03_features
```

No R, Rcpp, RcppArmadillo, or RcppParallel installation is required.

## Minimal use

```fortran
use rfast

real(dp) :: x(5)
type(mle_result) :: fit

x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
print *, mean_r(x), median_r(x), variance_r(x)
fit = normal_mle(x)
print *, fit%param
```

## Compatibility notes

R's dynamic types, data frames, formula/model-frame machinery, namespace manipulation, S3 behavior, NA/NaN policy, and nonconformable vector recycling are not emulated. Fortran callers use typed arrays and explicit matrices. Parallel Rcpp/RcppParallel switches are not reproduced; independent column operations can be parallelized by the caller.

The multinomial solver reports a nonzero status for singular/separated designs instead of hiding the rank failure. `mv_eeltest2` currently implements the upstream standard chi-square calibration; the optional James/F calibration variants remain pending. Upstream `omp` is complete computationally, but its rank-specific Fortran API is split into `omp_glm`, `omp_multivariate`, and `omp_multinomial`. `pc_skeleton` currently covers Pearson/Spearman analytic testing; categorical and permutation modes remain pending.

See `API_MAP.md` and `PORTING_NOTES.md` for exact coverage and limitations.
