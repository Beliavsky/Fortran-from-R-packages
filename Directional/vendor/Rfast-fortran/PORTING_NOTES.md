# Porting notes

## Scope

Rfast 2.1.5.2 is a library-scale package with 444 unique NAMESPACE exports. Version 0.3.0 ports 257 exports/consolidated computational interfaces, one partial export, and identifies every remaining export in `API_MAP.md`. It is an expanded computational-core port, not a claim of complete R-package parity.

## Array orientation

The Fortran API uses ordinary mathematical orientation: observations are rows and variables are columns. Several Rfast C++ routines receive transposed R matrices internally for Rcpp efficiency; those implementation transposes are not exposed in the Fortran API.

## Random numbers

`zigg-fortran` is a local FPM dependency. `set_seed`, `rnorm_vector`, `rexp_vector`, and `runif_vector` use its generator. Distribution-specific generators are built on that stream. This keeps the numerical library independent of the R runtime.

## Numerical linear algebra

The core is self-contained and does not require BLAS/LAPACK. Cholesky, Gaussian-elimination solves/inversion, and a Jacobi symmetric eigensolver are included. Large dense problems would benefit from an optional BLAS/LAPACK backend in a later performance release.

## v0.2 score tests

The ten `score.*regs` exports are represented by native Fortran routines returning statistics and p-values. They reuse the corresponding scalar null-model MLE kernels already in the port and evaluate candidate columns without reproducing R formula/model-frame machinery.

## Multivariate MLEs

`mvnorm_mle` and `mvlognormal_mle` use the exact sample MLE formulas. `mvt_mle` uses an iterative Student-t location/scatter update. `dirichlet_mle` and `inverse_dirichlet_mle` use Newton updates with digamma/trigamma support. Nonconvergence or singular updates are reported through the result status instead of raising an R condition.

## Empirical likelihood

`el_test1`/`el_test2` solve the standard empirical-likelihood constraints directly. `eel_test1`/`eel_test2` use exponential tilting, and the multivariate EEL routines solve the vector moment equations with Newton updates. `mv_eeltest2` currently implements the standard chi-square calibration (the upstream `R=0` path); optional James/F calibration variants are still pending.

## Repeated measures and variance components

`rm_anova`, `rm_anovas`, and `rm_lines` use native Fortran array shapes rather than R's flattened matrix conventions.

For balanced random-intercept data, `colvarcomps_mle` dispatches to the upstream closed-form balanced calculation. For unbalanced groups, it translates the Rfast C++ iteration based on `gold_rat3`, including the upstream variance-ratio search interval `[0,50]`. This is intentionally an algorithm-compatibility choice; it need not equal an unconstrained generic mixed-model optimizer on cases whose optimum lies outside that interval. `varcomps_mle` is the single-response convenience wrapper.

## Regression and selection

The v0.2 regression module adds Gamma and inverse-Gaussian log-link regressions, quasi-Poisson, proportion/logistic quasi-likelihood, baseline-category multinomial regression, and multivariate spatial-median regression. Perfect separation or singular multinomial information matrices return a nonzero status.

`ompr` implements OMP-style forward selection for the BIC, SSE, and adjusted-R2 criteria. `bic_corfsreg` implements partial-correlation BIC forward selection. `omp_glm` currently covers the logistic and Poisson upstream branches; quasi-Poisson, quasi-binomial, normal-log, Gamma, Weibull, and multivariate `omp` branches remain pending.


## v0.3 regressions and random intercepts

`normlog_regression` and `weibull_regression` translate the upstream Newton updates. Independent numerical objective checks are included in the v0.3 regression test. `tobit_mle` is the upstream left-censored-at-zero location/scale MLE, and `ordinal_mle` reproduces the empirical-threshold logit/probit estimator.

The `rint_*` routines translate Rfast's own random-intercept algorithm, including the bounded `gold_rat3` variance-ratio optimization. This deliberately targets upstream compatibility rather than replacing it with a generic mixed-model optimizer. The supplied group IDs may be unbalanced.

## v0.3 OMP and PC skeleton

All upstream `omp` families are now implemented. Because Fortran cannot overload one routine on both real vector, integer vector, and real matrix responses without a generic wrapper ambiguity in this API, the implementation is exposed as `omp_glm`, `omp_multivariate`, and `omp_multinomial`. The v0.3 code also fixes the v0.2 selection-loop fidelity issue by checking the first fitted model's improvement before selecting a second variable. For multivariate OMP, the mathematically consistent decrease in the covariance log-determinant criterion is used; the supplied R source tests the opposite sign in its loop condition.

`pc_skeleton` implements the Pearson and Spearman analytic (`R=1`) modes with Fisher-z conditional-correlation tests, separation sets, and arbitrary feasible conditioning order. Upstream categorical G2 mode and permutation-correlation mode (`R>1`) are still pending, so `pc.skel` remains marked partial. Missing-value imputation from the R wrapper is also not implicit; Fortran callers should preprocess missing data explicitly.

## Statistical conventions

- `variance_r` defaults to sample variance (`n-1` denominator); pass `population=.true.` for the population/MLE denominator.
- Normal and lognormal MLE variances use the population denominator.
- The two-sample `energy_distance` follows the exact scaling used by upstream `Rfast::edist`.
- `binomial_mle` implements the common fixed-number-of-trials branch; upstream's unknown-trials Newton branch remains pending.
- `vmf_mle` supplies mean direction and the Rfast-style kappa approximation; arbitrary-dimensional Bessel normalization/log-likelihood remains pending.

## R-only behavior intentionally omitted

Namespace/hash/iterator utilities, data-frame conversion, package checking/source helpers, S3 dispatch, plotting, and formula/model-frame construction are R infrastructure rather than standalone numerical algorithms and are not recreated.

## Advanced computational routines still pending

The largest remaining clusters include full forward/backward model-selection families beyond OMP, categorical/permutation PC-skeleton modes, high-dimensional outlier algorithms, specialized directional distributions, robust/specialized regression families, and several less-common MLE/test families. Each is listed individually in `API_MAP.md`.
