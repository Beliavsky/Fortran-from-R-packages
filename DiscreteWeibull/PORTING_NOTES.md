# Porting notes

## Scope

This is a complete translation of the computational API of
`DiscreteWeibull` 1.1. There is no plotting code in the package source;
graphics occur only in documentation examples.

## Rsolnp dependency

The supplied `Rsolnp-fortran` 2.0.1 translation is vendored and used by the
type-I method-of-moments (`method="M"`) fit, matching the upstream package's
use of `solnp`. Callback data are carried through `solnp_problem%data`, so the
Fortran implementation does not rely on module-global observations.

## Parameter domains

The implementation enforces the domains documented by the package:

- type I: `0 <= q < 1`, `beta > 0`;
- type III: `c > 0`, `beta >= -1`.

ML fits use unconstrained transformed coordinates internally and therefore
cannot step outside these domains.

## Numerical stability

PMFs are evaluated in log space with stable log-difference formulas. Type-I
quantiles explicitly handle probabilities 0 and 1. For type III at
`beta=-1`, quantile bracketing uses harmonic-number evaluation rather than the
upstream unbounded linear loop.

The finite-support/continuous-Weibull approximations in `Edweibull` and
`E2dweibull` are retained for compatibility with the upstream package.

## Upstream hazard typo

The supplied R source defines

`hdweibull3 <- function(x,c,beta) 1-exp(-c(x+1)^beta)`.

Because `c(x+1)` is parsed as a call to R's `c()` function, the parameter `c`
does not multiply the power as intended. The package documentation defines
the type-III cumulative hazard with parameter `c`, so the Fortran translation
uses the mathematically consistent hazard

`1 - exp(-c*(x+1)^beta)`.

## Fisher information

The R function differentiates the log mass symbolically and sums the observed
second derivatives. The Fortran implementation computes the same observed
information matrix as the central finite-difference Hessian of the average
negative log-likelihood at the ML estimate. The returned inverse matches the
upstream convention: divide it by sample size to obtain the asymptotic
covariance of the ML estimator.

## RNG

The inverse-transform algorithms are retained, but Fortran's intrinsic RNG is
used; seeded streams therefore do not match R's `runif()` bit-for-bit.
