# API overview

Import the facade with:

```fortran
use survey
```

The facade re-exports the public procedures and result/design types from the
modules below.

- `survey_design`: `make_design`, `validate_design`, `design_weights`, `design_degf`
- `survey_estimators`: totals, means, ratios, covariance, cross-products, CDF/table/kappa primitives
- `survey_replicates`: JK1/JKn/BRR/bootstrap constructors and replicate estimators
- `survey_calibration`: calibration, post-stratification, raking, trimming
- `survey_quantiles`: weighted/Taylor/replicate quantiles and qrule constants
- `survey_glm`: survey GLM, replicate GLM, prediction and pseudo-R2 helpers
- `survey_chisq`: Wald and Rao-Scott contingency tests
- `survey_ivreg`: survey and replicate 2SLS
- `survey_pps`: joint-inclusion approximations and HT/Yates-Grundy variance
- `survey_inference`: t tests, rank tests, proportion CIs, contrasts
- `survey_mle`: generic maximum pseudo-likelihood callbacks
- `survey_nls`: generic nonlinear survey regression callbacks
- `survey_survival`: KM, Cox, log-rank, and parametric AFT survey wrappers
- `survey_multivariate`: Cronbach alpha, weighted correlation, PCA
- `survey_special`: incomplete-beta/gamma and F/chi-square tail helpers

## Core design type

`survey_design_t` stores numeric weights plus stage-by-stage PSU, stratum, and
FPC arrays. `make_design` validates common shapes and creates the design.
Estimator routines return `svystat_t` or a specialized result type containing
estimates, covariance matrices, degrees of freedom, and convergence metadata
where applicable.

## Callback APIs

`survey_mle` and `survey_nls` use Fortran procedure callbacks rather than R
formulas. The caller owns parameter ordering and model-matrix construction.
This is deliberate: it exposes the numerical contract directly and avoids an
embedded expression parser.

`survey_mle` uses module-level callback dispatch while `minqa` is active, so a
single process should not run simultaneous `svy_mle` fits in separate threads.
