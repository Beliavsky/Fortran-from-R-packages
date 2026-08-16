# VGAM-fortran

A modern Fortran 2018 / FPM computational port of numerical functionality from
Thomas W. Yee's R package **VGAM 1.1-14**.

This is a numerical library, not an emulation of R's S3/S4 object system. It
keeps model fitting, link functions, distributions, special functions,
constraint matrices, reduced-rank fitting, categorical regression, count
models, time-series likelihoods, and spline-based additive fitting in explicit
typed Fortran APIs. Plotting and other graphics code are omitted.

The supplied `splines-fortran-v0.1.0` translation is vendored and used by the
additive/spline routines rather than reimplementing R's `splines` dependency.

## Implemented computational areas

- Stable link functions and inverse links: identity, log, logit, probit,
  complementary log-log, log-log, cauchit, reciprocal, square-root, Fisher-z,
  and negative variants, with derivatives.
- Portable special functions used by VGAM-style distributions: normal CDF and
  quantile, incomplete gamma/beta, polygamma, Hurwitz/Riemann zeta, Lambert W,
  exponential integrals, Lerch phi series, generalized harmonic numbers, and
  Kendall tau.
- Standard density/CDF/quantile/random-number routines for normal, logistic,
  exponential, gamma, beta, Poisson, binomial, and negative-binomial laws.
- Dirichlet density/simulation and joint Dirichlet regression with log-shape
  predictors, optional parallel slopes, fitted composition means, covariance, and
  exact family-specific expected-information matrices on shape/log-shape scales.
- VGAM distribution functions for Gumbel, Frechet, Rayleigh/generalized
  Rayleigh, Pareto I/II/IV, beta-binomial, Dirichlet-multinomial, zeta/Zipf,
  inverse Gaussian, exponential-geometric/Poisson/logarithmic, triangular,
  GEV, generalized Pareto, Laplace, and Kumaraswamy distributions.
- Actuarial distributions: Gompertz, Makeham, Perks, Lindley, Nakagami,
  Maxwell, Benini, Levy, and Skellam.
- GLM/VGLM-style IRLS for Gaussian, Poisson, binomial, Gamma, and inverse
  Gaussian responses, with weights, offsets, penalties, covariance, standard
  errors, deviance, log-likelihood, AIC, and prediction.
- Joint baseline-category multinomial logit fitting and cumulative-link ordinal
  regression.
- Beta regression, negative-binomial regression, and zero-inflated Poisson.
- Zero-truncated Poisson, hurdle Poisson, hurdle negative binomial, and
  zero-inflated negative-binomial regression with separate count/zero designs.
- VGAM-style linear coefficient constraints for multi-response independent
  vector GLMs, including identity and parallel-constraint helpers.
- Alternating reduced-rank VGLM fitting for independent supported families,
  with an unrestricted predictor block and factorization
  `B_reduced = C A^T`. The default leaves an intercept column unrestricted.
- Doubly constrained reduced-rank VGLM fitting with explicit per-latent
  `H.A` column-space constraints and per-reduced-predictor `H.C` constraints,
  including mixed supported families, weights, offsets, prediction, and rank
  diagnostics.
- Quadratic RR-VGLM fitting with the full symmetric latent quadratic form
  `eta_j = X1*b_j + A_j*z + z^T Q_j z`, cross-latent curvature, optional
  `Dzero`-style linear-only responses, prediction, and latent optimum/curvature
  extraction.
- CQO computational API built on QRR-VGLM: fitting, latent response surfaces,
  response-to-latent calibration, optimum/curvature extraction, and minimum-norm
  reconstruction of reduced environmental predictors from calibrated scores.
- Rank-1 constrained additive ordination (CAO), matching the principal upstream
  current restriction, with cubic penalized-spline response curves and alternating
  estimation of the canonical environmental coefficient vector.
- A reusable finite-support GAITD mass-transform kernel supporting truncation,
  fixed altered masses, additive inflated masses, and multiplicative deflation,
  with Poisson and negative-binomial base wrappers, CDFs, quantiles, RNG, and
  moments.
- Covariate-dependent GAITD Poisson and negative-binomial regression with
  multinomial-logit baseline/special-mass probabilities, altered and inflated
  points, fixed truncation, likelihood fitting, covariance, moments/prediction,
  and separate mean/mass design matrices.
- GAITD MLM semantics for Poisson and negative-binomial parents: exact direct
  altered probabilities, additive inflation, additive deflation, truncation,
  PMF/CDF/quantiles/moments, plus covariate-dependent MLM regression with
  altered/inflated/deflated special points.
- GAITD outer-distribution `a.mix`/`i.mix`/`d.mix` constructors for Poisson and
  negative-binomial parents. A total special mass is distributed across its named
  support points according to a separately parameterized outer Poisson/NB law,
  and mix components can coexist with direct MLM masses and truncation.
- Covariate-dependent GAITD outer-mix regression for Poisson and negative-binomial
  parents. Parent means, mix masses, and outer-distribution means have explicit
  design matrices. The original NB fit estimates a shared dispersion, while v0.8
  also provides separate log-dispersion predictors for the parent and each outer
  altered/inflated/deflated NB law.
- Bivariate copula kernels for Clayton, Frank, Farlie-Gumbel-Morgenstern (FGM),
  Gaussian, Plackett, and Ali-Mikhail-Haq (AMH) copulas, including VGAM-compatible
  density/CDF parameterizations, random generators, and one-parameter dependence
  regression.
- Student-t numerical layer: univariate density/CDF/quantile/RNG, the VGAM
  bivariate Student-t density and likelihood fit with covariate-dependent degrees
  of freedom and correlation, plus `dbistudenttcop`-compatible Student-t copula
  density, RNG, and correlation regression with fixed or estimated degrees of freedom.
- Zero/one-altered beta `d/p/q/r` functions and zero/one-inflated beta-binomial
  PMF/CDF/RNG kernels, preserving the endpoint-mass semantics of the upstream
  computational helpers.
- Positive-normal and positive-geometric distribution helpers; zero-altered
  Poisson, negative-binomial, geometric, and binomial `d/p/q/r` helpers; and
  zero-inflated/deflated binomial and geometric helpers with the upstream
  admissible deflation limits.
- Zero-altered Poisson, negative-binomial, geometric, and binomial regression with
  separate parent-parameter and observed-zero predictors; NB dispersion is jointly
  estimated.
- Censored normal, Poisson, exponential, and Rayleigh regression with exact,
  left-, right-, and interval-censored observations, numerical covariance, AIC,
  and array-oriented prediction.
- Tobit distribution/regression with exact endpoint masses, separate latent-mean
  and log-scale designs, plus generalized folded-normal d/p/q/r and regression.
- Positive/zero-truncated Poisson helpers and positive negative-binomial
  d/p/q/r/regression with jointly estimated NB dispersion.
- Zero/one-altered beta regression with separate beta-mean, precision, and zero/one
  mass predictor blocks constrained by a numerically stable softmax.
- Additional multivariate-family computations: bivariate normal density/RNG and
  regression, VGAM bivariate logistic density/CDF/RNG and regression, Freund
  (1961) bivariate exponential density/RNG/regression, trivariate normal
  density/RNG/MLE, and the VGAM FGM-exponential (`bifgmexp`) density/CDF/RNG and
  dependence regression.
- Information/constraint utilities for numerical observed information, score
  outer-product information, constrained information projection, and covariance
  lifting from free to full coefficient spaces.
- Exact Yeo-Johnson transform/inverse and lambda derivatives, plus a
  likelihood-based normal/Yeo-Johnson regression fit with quantile prediction.
- Three-predictor LMS/Yeo-Johnson likelihood fitting with separate design
  matrices for `lambda(x)`, transformed location `mu(x)`, and `log sigma(x)`,
  compatible with spline bases and with original-scale quantile prediction.
- Gaussian AR(1) exact/conditional likelihood fitting and forecasting.
- Nested reduced-rank autoregression (RRAR) with non-increasing lag ranks,
  a shared identified left subspace, lag-specific right factors, concentrated
  Gaussian likelihood, innovation covariance, parameter covariance, and recursive
  forecasting.
- GARMA(p,0) fitting for identity, log/reciprocal, logit, probit, cloglog, and
  cauchit links, with the upstream observation-driven AR predictor structure,
  numerical covariance, and recursive forecasting.
- Penalized B-spline / natural-spline additive fits for Gaussian, Poisson, and
  binomial responses, using the supplied `splines-fortran` backend.
- Dense numerical linear algebra and a self-contained BFGS/numerical-Hessian
  optimizer used by multi-parameter families.

## Deliberately not translated

Plotting/graphics, R formula parsing, model frames, S3/S4 dispatch, printing,
summary formatting, smart-prediction bookkeeping tied to R environments, and
R-specific methods whose main purpose is integration with other R classes.

## Important v0.9 scope note

VGAM is very large. Version 0.9 substantially extends the numerical port, but it
is **not yet a one-for-one implementation of every advanced VGAM family**.
DRR now accepts explicit `H.A`/`H.C` linear subspace constraints, CQO has fitting
and calibration utilities, rank-1 CAO has a spline-based alternating fitter, RRAR
is represented as a nested reduced-rank VAR, and GAITD now has ordinary special-mass
regression, direct MLM semantics, and explicit outer-distribution mix constructors.
Six one-parameter copula families plus Student-t copula computations are
also included. Version 0.7 additionally includes direct outer-mix regression,
four censored-regression families, and bivariate normal/logistic/Freund models.
Version 0.8 adds separate outer-NB dispersion regression, common positive and
zero-altered count distributions/regressions, trivariate normal fitting, and
FGM-exponential dependence modeling. These are array-oriented numerical APIs,
not clones of the R control/formula/S4 layers.
Version 0.9 adds joint Dirichlet regression with an analytic family-specific EIM,
Tobit and generalized folded-normal likelihood models, positive NB regression,
and zero/one-altered beta regression. These additions remain explicit numerical
APIs rather than R family-object clones.

Remaining large areas include the complete upstream DRR corner/summary covariance
bookkeeping, higher-level CQO/CAO control and calibration methods, the full
generated GAITD system (especially exact generated EIM/control code for every
wrapper and the complete wrapper catalog), the remaining bivariate/copula
catalog, and remaining specialized altered/inflated/censored families and R family-object wrappers. See
`docs/TRANSLATION_NOTES.md`.

## Build

With FPM:

```text
fpm build
fpm test
fpm run --example
```

The sources use standard free-form Fortran 2018 and do not require compiler
extensions or nonstandard line-length flags.

## Minimal use

```fortran
program demo
   use vgam
   implicit none
   real(dp) :: x(5,2), y(5)
   type(vglm_result_t) :: fit
   integer :: i

   do i=1,5
      x(i,:) = [1.0_dp, real(i-3,dp)]
   end do
   y = [0.2_dp, 0.4_dp, 0.5_dp, 0.8_dp, 0.9_dp]

   call fit_binomial(y,x,fit)
   print *, fit%coefficients
end program demo
```

See `example/rrvglm_v02.f90`, `example/qrr_garma_v03.f90`,
`example/drr_cao_rrar_v04.f90`, `example/gaitd_copula_v05.f90`, and
`example/student_mix_v06.f90`, `example/censored_bivariate_v07.f90`, and
`example/zero_altered_multivariate_v08.f90` for reduced-rank, ordination,
time-series, GAITD, copula, Student-t, censored, altered-count, and multivariate
APIs.

## License and attribution

The original VGAM package declares GPL-3 in `DESCRIPTION`; its included
`LICENCE.note` contains additional provenance and licensing information. Those
files are preserved verbatim under `original/` and the top-level `LICENSE`
contains GPL version 3.

The vendored spline translation is GPL-2.0-or-later and retains its own license
and attribution under `vendor/splines/`. Combining it into this derivative
package under GPL-3.0 is license-compatible.

No endorsement by the original VGAM authors, R Core, or the authors of the
splines package is implied.
