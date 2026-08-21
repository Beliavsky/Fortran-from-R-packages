# Translation coverage

## Exported upstream APIs with direct numerical counterparts

| R API | Fortran counterpart |
|---|---|
| `boot` | `bootstrap_weighted` |
| `boot.array`, `freq.array` | resampling arrays + `frequency_array` |
| `norm.ci` | `normal_ci` |
| basic bootstrap CI | `basic_ci` |
| studentized CI | `studentized_ci` |
| percentile CI | `percentile_ci` |
| BCa CI | `bca_ci` |
| `abc.ci` | `abc_ci` |
| `empinf` | infinitesimal/delete-one/positive/regression influence routines |
| `linear.approx` | `linear_approximation` |
| `exp.tilt` | `exponential_tilt` |
| `imp.weights` | `importance_weights` |
| `imp.moments` | `importance_moments` |
| `imp.quantile` | `importance_quantile` |
| `imp.prob` | `importance_probability` |
| `corr` | `weighted_corr` |
| `var.linear` | `var_linear` |
| `k3.linear` | `k3_linear` |
| `cum3` | `cum3` |
| `logit`, `inv.logit` | `logit`, `inv_logit` |
| `simplex` | `simplex_solve` |
| `saddle` (simple multinomial) | `multinomial_saddlepoint` |
| `envelope` | `confidence_envelope` |
| `smooth.f` | `smooth_frequencies` |
| `EL.profile` numerical kernel | `empirical_loglikelihood` |
| `EEF.profile` numerical kernel | `eef_loglikelihood` |
| `lik.CI` | `likelihood_ci` |
| `nested.corr` | `nested_correlation` |
| `glm.diag` numerical formulas | `glm_diagnostics`, scale helpers |
| `tsboot` fixed/geometric resampling core | block-index generators |
| `censboot` case/product-limit core | `boot_censored` primitives |

## Higher-level R wrappers represented by lower-level Fortran primitives

- `boot.ci`
- `tilt.boot`
- `tsboot`
- `control` (partially; see below)

## Remaining numerical targets

These are not graphics, but still depend strongly on R model/spline objects or
specialized interfaces and are not complete in v0.1.0:

- `cv.glm` model refitting/prediction wrapper
- full `censboot` model/conditional/weird orchestration, especially `coxph`
  model-based failure-time sampling
- conditional Poisson and binary modes of `saddle`
- `saddle.distn` smoothing-spline distribution construction
- the full `control` quantile path that consumes a `saddle.distn` spline
- the highest-level replay/reconstruction semantics of R `boot.array`

## Intentionally omitted presentation code

- `plot.boot`
- `print.boot`
- `print.bootci`
- `print.simplex`
- `print.saddle.distn`
- `lines.saddle.distn`
- `glm.diag.plots`
- graphical part of `jack.after.boot`
