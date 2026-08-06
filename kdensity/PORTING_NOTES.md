# Porting notes

## Numerical architecture

The R package returns a callable closure with an attached environment. The
Fortran port returns `type(kdensity_fit)`, which stores the sample, fitted start,
kernel, bandwidth, support, and normalization constant. Type-bound procedures
perform scalar and vector density evaluation.

Adaptive Simpson integration is used for normalization. Finite intervals use a
sine-squared endpoint transformation, which handles integrable beta-kernel and
beta-start endpoint singularities. Infinite intervals use rational or tangent
transformations.

## Parametric starts

The R package discovers starts dynamically from `univariateML`. That dependency
was not supplied with this translation. The port independently implements the
common starts used by the package documentation and tests and exposes the same
extensibility through `type(kd_start)` procedure pointers.

Normal, lognormal, exponential, Weibull, inverse-Gaussian, Laplace, and Pareto
starts use direct likelihood or closed-form estimators. Several remaining starts
use stable moment/quantile estimators. Consequently, fitted start parameters may
differ slightly from `univariateML` maximum-likelihood estimates.

## Bandwidth selectors

The JH, RHE, HS, nrd0, and nrd formulas are direct translations. Generic UCV
uses leave-one-out start refits, normalized squared-density quadrature, and a
golden-section search.

R's `bw.bcv` and `bw.SJ` are large internal algorithms not implemented in the
upstream `kdensity` source. The compatibility names use documented reference
approximations. Users needing exact R parity should supply a numerical bandwidth
explicitly.

## Source behavior and corrections

- The upstream object stores `sum(density(x))` in a field called `logLik`.
  The Fortran port stores the mathematically correct `sum(log(density(x)))`.
- The upstream default with no explicit kernel or support is Gaussian even when
  a bounded parametric start is selected. The Fortran default preserves this.
- `bw = Inf` is represented by a sufficiently large IEEE finite value in the
  typed options and returns the fitted parametric start directly.
- Endpoint transformations prevent false normalization explosions for
  integrable singular densities.

## Custom procedures

A custom start supplies a pure scalar density and an estimator subroutine. A
custom kernel supplies a pure scalar `(y, x, h)` function. Assign those procedure
pointers to `type(kd_start)` and `type(kd_kernel)`, then call
`fit_kdensity_custom`.
