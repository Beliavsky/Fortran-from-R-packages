# Porting notes

## Interface strategy

R's `boot()` accepts arbitrary R objects and functions and stores replayable
calls, RNG state, formulas, and S3 metadata.  The Fortran port exposes the
numerical operations directly.  Statistics are supplied through explicit
procedure callbacks and data are numerical arrays.

`bootstrap_weighted` currently targets the important nonparametric weighted
statistic path.  It returns the observed statistic, replicate statistics,
resampling indices, and frequency arrays.

## Numerical compatibility

- `norm_inter` follows the upstream interpolation rule on the normal-quantile
  scale rather than using an ordinary type-7 sample quantile.
- `imp.moments` preserves the upstream special case in which constant
  importance weights produce the ordinary sample variance.
- `linear_approximation` uses the same stratum-wise `frequency * influence / n_s`
  expression as R `linear.approx`.
- `exponential_tilt` solves the same tilted-mean equation, using Newton updates
  rather than R `optim(BFGS)` on the squared mismatch.
- the LP routine uses a revised two-phase simplex implementation rather than the
  upstream tableau implementation.  It solves the same nonnegative-variable
  problem with `<=`, `>=`, and equality constraints.
- the simple multinomial saddlepoint implementation solves the cumulant
  saddlepoint equations by Newton iteration and implements the same
  Barndorff-Nielsen/Lugannani-Rice one-dimensional CDF corrections.
- Fortran `random_number` is used for stochastic routines.  Statistical behavior
  is preserved, but random streams are not expected to reproduce R's RNG bit for
  bit.

## Deliberate API differences

- plotting and S3 print/plot methods are omitted
- `boot.ci` is represented by its numerical component routines instead of one
  R-object dispatcher
- `tilt.boot` is represented by the bootstrap, influence, exponential-tilt,
  importance-resampling, and smoothing primitives
- `tsboot` is represented by its fixed/geometric block index generators; model
  simulation remains a caller-supplied generator task
- `censboot` currently exposes case resampling plus product-limit and conditional
  censoring primitives.  The full R `survfit`/`coxph` object wrapper and Cox
  model-based resampling are not yet implemented
- `glm.diag` formulas are translated, but `cv.glm` is not: it depends on a
  general reusable GLM fit/predict object interface that is outside this package
- `saddle` currently implements the simple multinomial case.  Conditional
  Poisson/binary saddlepoints and the smoothing-spline `saddle.distn` wrapper
  remain targets for a later version
- `control` remains a higher-level composition target because its quantile path
  depends on `saddle.distn`

## Upstream source

The original package source is retained verbatim under `upstream/boot-1.3-32`
for license, attribution, and algorithm auditing.
