# Changelog

## 0.9.0 - 2026-08-14

- Added joint Dirichlet regression with log-shape predictors, optional parallel
  slopes, fitted shape/mean matrices, covariance, simulation, and exact analytic
  expected-information kernels for shape and log-shape parameterizations.
- Added VGAM-compatible Tobit d/p/q/r computations with lower/upper endpoint
  point masses and Tobit regression with separate latent-mean and log-scale
  design matrices.
- Added generalized folded-normal d/p/q/r computations and likelihood regression
  for separate mean and scale predictors and arbitrary positive `a1`/`a2` folds.
- Added positive/zero-truncated Poisson d/p/q/r helpers and positive
  negative-binomial d/p/q/r plus regression with jointly estimated dispersion.
- Added zero/one-altered beta regression with separate beta-mean, precision, and
  endpoint-mass predictor blocks using a stable three-category softmax.
- Added `test_v09` and `example/dirichlet_tobit_v09.f90`; all v0.1-v0.9 tests
  and all examples pass under Fortran runtime checking.

## 0.8.0 - 2026-08-14

- Added `fit_gaitd_mix_nb_dispersion_regression`, allowing separate log-dispersion
  predictors for the negative-binomial parent and each of the altered, inflated,
  and deflated outer NB distributions while retaining the exact restricted-support
  GAITD mix normalization.
- Added positive-normal and positive-geometric distribution helpers and direct
  zero-altered Poisson, negative-binomial, geometric, and binomial d/p/q/r APIs.
  Added zero-inflated/deflated binomial and geometric d/p/q/r helpers with the
  upstream admissible negative-mass limits.
- Added zero-altered Poisson, negative-binomial, geometric, and binomial regression
  with separate parent and observed-zero predictors; the NB model also estimates
  dispersion.
- Added trivariate-normal density/RNG/MLE with unrestricted means, scales, and
  three correlations, including positive-definiteness checks and covariance output.
- Added the VGAM `bifgmexp` FGM-exponential density/CDF/RNG and dependence regression,
  plus a Kendall-tau computational helper.
- Added `test_v08` and `example/zero_altered_multivariate_v08.f90`. A clean build
  passes all v0.1-v0.8 tests and all examples under Fortran runtime checking.

## 0.7.0 - 2026-08-13

- Added covariate-dependent GAITD outer-mix regression for Poisson and
  negative-binomial parents. Parent mean, special-mass logits, and the three
  outer-distribution means have separate design matrices; the NB implementation
  uses one estimated size shared by parent and outer NB components.
- Added likelihood-based censored regression for normal, Poisson, exponential,
  and Rayleigh families with exact, left-, right-, and interval-censored
  observations, prediction, AIC, and numerical-Hessian covariance.
- Added bivariate normal density/RNG/regression with covariate-dependent means,
  scales, and correlation; VGAM bivariate logistic density/CDF/RNG/regression;
  and Freund (1961) bivariate exponential density/RNG/regression.
- Added information/constraint helpers for observed numerical information,
  score outer-product information, constrained information projection, and
  lifting free-parameter covariance back to the full coefficient space.
- Replaced an unstable exponential cap in the bivariate-logistic likelihood
  with a log-sum-exp calculation after a boundary-regression test exposed the
  numerical failure mode.
- Added `test_v07` and `example/censored_bivariate_v07.f90`; all previous tests
  and examples continue to pass with runtime checking.

## 0.6.0 - 2026-08-13

- Added GAITD outer-distribution `a.mix`/`i.mix`/`d.mix` constructors for Poisson
  and negative-binomial parents, including coexistence with direct MLM masses
  and truncation and exact restricted outer-distribution point weighting.
- Added Ali-Mikhail-Haq copula density/CDF/RNG and one-parameter regression.
- Added self-contained Student-t density/CDF/quantile/RNG, bivariate Student-t
  density/simulation/likelihood fitting, and Student-t copula density/RNG/fitting.
- Added zero/one-altered beta d/p/q/r and zero/one-inflated beta-binomial
  probability/CDF/RNG computational helpers.
- Added `test_v06` and `example/student_mix_v06.f90`; all prior tests/examples
  remain passing under Fortran runtime checking.

## 0.5.0 - 2026-08-13

- Added direct GAITD MLM Poisson/negative-binomial constructors with upstream-style
  direct altered masses, additive inflation, additive deflation, truncation,
  normalization, CDF/quantiles/RNG, and moments.
- Added `fit_gaitd_mlm_poisson_regression` and `fit_gaitd_mlm_nb_regression` with
  covariate-dependent direct special probabilities and pointwise validity checks.
- Added Clayton, Frank, FGM, Gaussian, and Plackett copula density/CDF kernels and
  random generators, preserving VGAM parameter conventions.
- Added `fit_copula_regression` for covariate-dependent one-parameter dependence
  models, covariance/AIC output, and dependence-aware initialization.
- Added `test_v05` and `example/gaitd_copula_v05.f90`; all earlier tests remain
  passing under Fortran runtime checking.

## 0.4.0 - 2026-08-13

- Added `fit_drrvglm` with explicit per-latent `H.A` loading constraints and
  per-reduced-predictor `H.C` latent-coefficient constraints. The constrained
  subspaces are optimized directly rather than imposed after an RR fit.
- Added nested reduced-rank autoregression (`fit_rrar`) with non-increasing lag
  ranks, a shared identified left subspace, lag-specific right factors,
  concentrated Gaussian likelihood, innovation covariance, optional parameter
  covariance, transformed series for the full-leading-rank case, and forecasts.
- Added a CQO API around the full QRR core: fitting, latent response surfaces,
  reverse calibration from multivariate responses to latent scores, and
  minimum-norm reconstruction of the reduced environmental predictor vector.
- Added rank-1 CAO with cubic penalized-spline response curves and alternating
  estimation of the canonical environmental direction, using the vendored
  `splines-fortran` backend.
- Added GAITD Poisson and negative-binomial regression with separate mean and
  special-mass designs, multinomial-logit baseline/special probabilities,
  altered/inflated points, fixed truncation, covariance, and prediction.
- Added `test_v04` and `example/drr_cao_rrar_v04.f90`; all earlier tests remain
  passing under Fortran runtime checking.

## 0.3.0 - 2026-08-13

- Added full symmetric quadratic reduced-rank VGLM fitting. Each response uses
  `eta = X1*b + A*z + z^T Q z`, including cross-latent quadratic terms,
  `Dzero`-style responses with zero quadratic curvature, prediction, and
  latent-space optimum/curvature extraction.
- Added a three-predictor Yeo-Johnson LMS likelihood with separate design
  matrices for transformation (`lambda`), transformed location (`mu`), and
  log-scale (`sigma`), numerical covariance, and original-scale quantiles.
- Added GARMA(p,0), matching the upstream package's currently implemented
  zero-MA case, for identity, log/reciprocal, and common binary links, with
  numerical covariance and recursive forecasting.
- Added `test_v03` including rank-1 and rank-2 QRR tests; the rank-2 case
  explicitly exercises cross-quadratic latent terms.
- Added `example/qrr_garma_v03.f90` and updated the aggregate `use vgam` API.

## 0.2.0 - 2026-08-13

- Added an alternating reduced-rank VGLM core with unrestricted and reduced
  predictor blocks, latent scores/loadings, mixed supported families,
  offsets/weights, rank diagnostics, and prediction.
- Added zero-truncated Poisson, hurdle Poisson, hurdle negative-binomial, and
  zero-inflated negative-binomial density/fitting APIs.
- Added a general finite-support GAITD mass-transformation engine with
  truncation, altered masses, inflation, deflation, moments, CDF/quantile/RNG,
  and Poisson/negative-binomial wrappers.
- Added the upstream Yeo-Johnson transform/inverse and lambda derivatives, plus
  a likelihood-based normal/Yeo-Johnson regression model with quantile
  prediction.
- Added exact/conditional Gaussian AR(1) likelihood fitting and forecasting.
- Added `test_v02` and a reduced-rank example; retained passing v0.1 tests.

## 0.1.0 - 2026-08-13

- Initial modern Fortran/FPM computational port of VGAM 1.1-14.
- Added links, special functions, random generators, distribution families,
  actuarial/extreme-value helpers, IRLS VGLM core, categorical/ordinal models,
  beta/negative-binomial/zero-inflated-Poisson regression, coefficient
  constraints, and penalized spline additive fitting.
- Vendored the user-supplied `splines-fortran-v0.1.0` dependency with its
  original translation license/notice retained.
- Added numerical and model integration tests plus a basic example.
- Plotting and R object/formula infrastructure intentionally omitted.
