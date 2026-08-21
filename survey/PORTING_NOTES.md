# Porting notes

## Scope

This is a modern Fortran/FPM translation of the computational core of R
`survey` 4.5. The port emphasizes statistical algorithms and explicit numeric
interfaces. R-specific object dispatch, formula parsing, model-frame creation,
DBI persistence, plotting, and presentation code are intentionally excluded.

This release is a substantial computational-core port, not a claim of complete
one-for-one coverage of every exported R function. Remaining numerical parity
targets are listed in `docs/TRANSLATION_COVERAGE.md`.

## Survey design representation

R's `survey.design` objects are represented by `survey_design_t`. Stages are
columns of integer PSU and stratum arrays. Sampling and population PSU counts
represent finite-population corrections explicitly. The recursive variance
engine follows the upstream multistage/onestage decomposition and includes the
ported lonely-PSU behavior.

The Fortran API expects callers to construct design/model matrices themselves.
This is the main semantic boundary between the R and Fortran interfaces.

## Replicate designs

`rep_design_t` stores full replicate weights. Constructors are provided for
JK1, JKn, BRR/Fay, and bootstrap replicates, together with replicate variance
and estimator routines. R-only automatic conversion, compressed-weight
storage, and several specialized replicate schemes are not yet complete.

## Calibration

The calibration solver translates the upstream generalized-raking Newton
iteration for linear/GREG, raking, bounded logit, and sinh calibration. Numeric
population totals/margins are passed directly rather than recovered from R
formula terms and factor levels.

## GLM and model fitting

The current `svy_glm` implements the most common canonical paths:
Gaussian/identity, binomial/logit, and Poisson/log. It uses survey estimating
functions and the design variance engine to construct a sandwich covariance;
it does not interpret survey weights as ordinary frequency weights.

`svy_mle` is a generic maximum pseudo-likelihood interface. It uses the supplied
`minqa` port for derivative-free optimization and the supplied `numDeriv` port
for numerical Hessians. When an observation-score callback is supplied, the
reported covariance is the survey sandwich covariance.

Because the supplied `minqa` callback API is scalar/global in style, `svy_mle`
uses module-level active callback state while an optimizer call is in progress.
Consequently simultaneous threaded `svy_mle` calls in one process are not
supported in v0.1.0.

`svy_nls` similarly exposes the nonlinear observation model as a procedure
callback and uses weighted Gauss-Newton plus the survey sandwich covariance.

## Survival integration

The user supplied a Fortran translation of R `survival`. Its top-level facade
was added to the active vendor copy so the survey wrappers can use its public
Cox/AFT routines. The user supplied current CRAN metadata identifying
`survival` 3.8-9 as LGPL (>= 2); the vendored metadata therefore records
`LGPL-2.0-or-later` and preserves the upstream author attribution.

The survey Cox wrapper fits the weighted Cox model through the supplied
survival library and computes design-robust covariance from observation-level
score contributions. Parametric survey survival regression uses the supplied
AFT likelihood and a numerical observation-score sandwich.

## Sparse Matrix dependencies

The supplied MatrixExtra translation is preserved under `vendor/` but is not
linked by the default FPM graph. The current implementation uses dense arrays.
This is deliberate both to keep the first parity release portable and to avoid
bringing the reference port's GPL-3-only linked choice into a default build
that also uses the supplied GPL-2-only `minqa` translation.

The separately supplied splines port is also reference-only because the active
survival translation already contains modules with the same `splines*` names.

## Numerical validation

`./scripts/strict_test.sh` compiles translated survey sources and the active
survival/numDeriv dependencies with strict GNU Fortran 2018 warning promotion
and runtime checking. The vendored minqa source is compiled without warning
promotion because its Powell-style implementation intentionally triggers
warnings for exact floating-point sentinel comparisons.

The current validation suite contains independent reference values for key
means/variances, GLM coefficients, Rao-Scott/Wald statistics, IV regression,
PPS variance, special-function probabilities, maximum pseudo-likelihood, and
nonlinear regression, plus integration/convergence tests for the supplied
survival translation.

## FPM availability

FPM was not installed in the build environment used for this port. The project
layout and `fpm.toml` were checked for TOML validity, and the exact source graph
was compiled and linked directly with GNU Fortran. On a machine with FPM, use
`fpm test`.

## v0.2.0 numerical parity expansion

Version 0.2.0 closes most of the high-value numerical gaps identified after the
initial computational-core port:

- `svyolr` is represented by an array/design-matrix cumulative-link fitter for
  logistic, probit, complementary-log-log/Gumbel, and cauchit links. Cutpoints
  are optimized through the same ordered reparameterization used upstream.
- `svyloglin` now has survey cell-probability covariance, loglinear fitting,
  and nested deviance/score tests based on the misspecification eigenvalues.
- `svyfactanal` now has ML factor fitting, effective-sample-size choices,
  Bartlett testing, and varimax rotation.
- native weighted chi-square/F Satterthwaite and saddlepoint tail calculations
  support Rao-Scott score/LRT inference. The optional upstream
  `CompQuadForm` numerical-integration path remains external.
- `Dcheck_strat`, multistage Dcheck construction, subset-aware Dcheck
  construction, and the two-phase Dcheck composition rule are translated.
  `multiphase_variance` directly implements the upstream sum of phase-specific
  HT covariance matrices and optionally residualizes influence functions in a
  calibration space before the HT calculation.
- two-frame constant and expected overlap estimators, totals/means, and
  independent-frame HT covariance are available through `survey_multiframe`.
- the Preston rescaled multistage bootstrap is available through `make_mrb`.
- pseudo-score, working Rao-Scott score, Wald term, and misspecified-likelihood
  ratio numerical tests are available without R formula/model-nesting code.

The remaining `svysmooth` algorithms intentionally are not replaced with a
look-alike smoother: upstream delegates local-polynomial smoothing and automatic
bandwidth selection to `KernSmooth`, and quantile smoothing to `quantreg`.
Likewise, `svysmoothArea` and `svysmoothUnit` currently delegate to the optional
`SUMMER` package. These are dependency translations rather than missing native
`survey` algorithms.
