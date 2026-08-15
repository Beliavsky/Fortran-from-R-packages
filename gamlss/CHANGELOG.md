# Changelog

## 0.9.0

- Added `marginal_predict_eta`, the matrix-first translation of upstream
  `getMarginal()` random-effect marginal prediction.
- Added the four upstream-compatible methods: `MARGINAL_INTEGRATE`,
  `MARGINAL_QFUNCTION`, `MARGINAL_RANDOM`, and `MARGINAL_NONE`.
- `MARGINAL_QFUNCTION` uses the exact 999-point `qnorm(0.001:0.999)` grid from
  upstream; `MARGINAL_RANDOM` defaults to 10,000 Gaussian draws per prediction.
- Added `get_marginal_random_intercept`, which removes each observation's
  fitted group effect before marginalizing the inverse link.
- Extended `random_intercept_result_t` with random-effect EDF and upstream-style
  `sigma_b = sqrt(sum(b**2)/edf)`.
- Added log-link analytic-oracle tests, identity-link tests, parameter-specific
  marginalization tests, object-adapter tests, and fitted random-scale tests.
- Added `example/v09_extended.f90` and retained upstream `R/random.R` under
  `upstream/reference-R-v09/` for provenance.


## 0.8.0

- Added `fit_gamlss_joint_random_effects_ais`, an adaptive deterministic
  quasi-Monte-Carlo importance-sampling marginal-likelihood path for joint
  cross-parameter random effects.
- Added group-specific posterior-mode/Hessian Gaussian proposals, Halton
  antithetic normal points, posterior covariance matrices, and effective sample
  size diagnostics.
- Supports combined random-effect dimensions up to 16 without tensor
  `order**dimension` growth.
- Added optional `refine_parameters=.true.` outer optimization of fixed effects
  and the full Cholesky-parameterized random-effect covariance under the AIS
  marginal likelihood.
- Added a damped-Newton inner mode solver so full AIS refinement does not nest
  the package BFGS optimizer recursively.
- Added `test_v08.f90`: a closed-form Gaussian random-intercept marginal
  likelihood oracle, a six-dimensional joint mu/sigma random-effect test, and
  an outer-refinement smoke test.
- Added `example/v08_extended.f90`.

## 0.7.0

- Added `mvn_rectangle_probability`, a deterministic Genz/Halton numerical
  evaluator for multivariate-normal rectangle probabilities, plus conditional
  Gaussian and MVN log-density helpers.
- Added `fit_gamlss_gaussian_copula_mixed` for discrete and mixed atomic/
  continuous GAMLSS margins.  Atomic observations use `F(y-)`/`F(y)` latent
  rectangles; continuous observations use exact conditional Gaussian density
  contributions.
- Added `family_cdf_left` and `family_observation_is_atom`, including endpoint
  atom handling for `BEINF`.
- Added `fit_gamlss_joint_random_effects_ghq`, a full marginal-likelihood
  Gauss-Hermite fitter for joint random-effect dimensions up to four, with
  full cross-parameter covariance and quadrature posterior effect means.
- Added `test_v07.f90` with a closed-form bivariate-normal rectangle oracle,
  grouped NBI copula recovery, BEINF atom checks, and a joint `mu`/`sigma` GHQ
  random-effects fit.
- Added `example/v07_extended.f90`.

## 0.6.0

- Added `fit_gamlss_gaussian_copula`, an exact Gaussian-copula joint likelihood
  for continuous GAMLSS margins with jointly optimized marginal regression and
  `nlme` correlation parameters.
- Added grouped temporal/spatial copula blocks, Gaussianized marginal scores,
  joint/marginal/copula log likelihood reporting and joint parameter covariance.
- Added `fit_gamlss_joint_random_effects`, estimating one full covariance over
  active distribution-parameter/random-term combinations, including
  cross-parameter covariance such as `Cov(b_mu,b_sigma)` and cross-slope terms.
- Added `test_v06` and `example/v06_extended.f90`.
- Retained the v0.5 correlated-RS estimator as a faster estimating-equation
  alternative to the new continuous-margin copula likelihood.


## 0.5.0

- Added `fit_gamlss_correlated_rs`, embedding `nlme` correlation matrices in
  non-Gaussian RS/Fisher working-response updates while retaining IRLS weights.
- Added `fit_gamlss_multi_random_effects` for simultaneous arbitrary grouped
  random-effect designs on multiple distribution parameters, with separate
  full within-parameter covariance matrices.
- Added full-backend `randomized_quantile_residuals_all`.
- Added `cross_validate_gamlss` with explicit K-fold held-out log scores and
  out-of-sample fitted parameters.
- Added `test_v05` and `example/v05_extended.f90`.


## 0.4.0

- Added `fit_gamlss_no_gls`, an exact Gaussian/NO adapter to the supplied `nlme`
  GLS backend with its translated correlation and variance-function structures.
- Added `fit_gamlss_multi_random_intercept` for simultaneous independent group
  random intercepts on any subset of `mu`, `sigma`, `nu`, and `tau`.
- Added `stepwise_gaic_parameter` for forward/backward/bidirectional GAIC/BIC
  selection on any fitted distribution parameter.
- Added numerical worm-plot coefficients, leverage-adjusted/Cook-style
  influence measures, and Jarque-Bera residual diagnostics.
- Refactored nested optimizer/integration callbacks in GAMLSS censoring and the
  vendored distribution backend so GNU Fortran no longer requires executable
  stack trampolines in validation builds.
- Added `test_v04` and `example/v04_extended.f90`.

## 0.3.0

- Added `fit_gamlss_random_effects` for arbitrary grouped random-intercept/slope
  designs with diagonal or unstructured correlated covariance updates.
- Added optional delayed-entry (`entry`) conditioning to generic censored
  GAMLSS likelihoods.
- Added `surv_interval2` and `(start,stop,event)` counting-process adapters.
- Added identifiable additive P-spline design/penalty composition.
- Added anisotropic 2D tensor-product P-splines.
- Added forward/backward/bidirectional matrix GAIC/BIC selection.
- Added case-resampling GAMLSS bootstrap coefficient/deviance sampling and
  percentile intervals.
- Added likelihood-ratio confidence intervals from coefficient profiles.
- Added `test_v03` and `example/v03_extended.f90`.
- Retained selected upstream random-effect, stepwise and profile R sources for
  v0.3 provenance.


## 0.2.0

- Added in-iteration ML-style scalar penalty/variance-ratio estimation to RS
  parameter updates.
- Added `fit_gamlss_random_intercept`, including optional Gaussian `nlme`
  variance-ratio initialization.
- Added generic 62-family CDF dispatch and exact/left/right/interval censored
  GAMLSS likelihood fitting.
- Added `surv_right_censoring` adapter for `(time,event)` data.
- Added fractional-polynomial basis/search/prediction.
- Added 1D/2D local-polynomial LOESS fitting and prediction.
- Added varying-coefficient P-spline construction.
- Added monotone P-spline fitting with active inequality penalties.
- Added iterative Lp `pcat` category fusion and fused-group extraction.
- Added matrix-based forward GAIC/BIC selection and coefficient profile
  likelihood.
- Added `test_v02` and `example/v02_extended.f90`.
- Retained selected upstream R sources used for v0.2 algorithm comparison.

## 0.1.0

Initial computational translation of `gamlss` 5.5-0.

- RS, CG and mixed fitting engines.
- Four parameter-specific design matrices, offsets, weights, fixes and steps.
- Quadratic penalties, EDF, AIC/GAIC/SBC and prediction.
- P-spline, natural-spline, cyclic, ridge and random-intercept helpers.
- Translation of compiled `genD.c` pairwise-difference kernel.
- Quantile residuals, residual summaries/ACF and model-comparison helpers.
- LMS/centile computation using BCCG/BCT/BCPE.
- Vendored supplied splines, nlme and survival translations.
- Vendored `gamlss.dist-fortran v0.3.0` distribution backend.
