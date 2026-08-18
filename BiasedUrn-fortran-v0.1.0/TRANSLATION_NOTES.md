# Translation notes

## Upstream

- Package: BiasedUrn
- Version: 2.0.12
- Date: 2024-06-16
- Author: Agner Fog
- License: GPL-3

The translation was made from the supplied `BiasedUrn-master` source tree.
The original R interface consists primarily of wrappers around Agner Fog's
C++ stochastic-library implementation.

## Coverage

All functions exported by the upstream `NAMESPACE` have a Fortran
computational counterpart:

- univariate Fisher density, CDF, quantile, random generation, mean,
  variance, mode, inverse odds, and inverse urn composition;
- univariate Wallenius equivalents;
- multivariate Fisher density, random generation, moments, inverse odds,
  and inverse urn composition;
- multivariate Wallenius equivalents;
- univariate and multivariate support limits.

R-specific vector recycling, SEXP allocation, R warnings/errors, dynamic
registration, and data-frame conversion are intentionally not translated.
Fortran callers use scalar functions or explicit arrays instead.

## Algorithm mapping

### Fisher distribution

The upstream C++ code has several table, inversion, and ratio-of-uniforms
paths selected for speed. The Fortran translation computes the same finite
support distribution with log-space normalization. For the multivariate
case, the normalizing coefficient is obtained by dynamic programming. This
is exact up to floating-point rounding and avoids explicit enumeration for
the density and RNG.

### Wallenius distribution

The upstream C++ implementation dynamically selects recursion, binomial
expansion, Laplace approximation, or numerical integration. The Fortran
translation uses the mathematical defining integral directly. It transforms
`t=exp(-u)`, locates the mode of the log integrand, truncates only negligible
tails according to the requested `precision`, and evaluates the scaled
integral with 32-point Gauss-Legendre quadrature on each side of the mode.

This is a deliberate algorithmic translation rather than a line-by-line
port of the upstream performance heuristics. Representative probabilities
were checked against independent implementations and agree to about
1e-11 absolute error or better at the default precision in the test region.

### Random generation

- Wallenius: exact sequential sampling without replacement with probabilities
  proportional to `remaining_i * odds_i`.
- Fisher: exact finite-support conditional sampling from log-space dynamic
  programming coefficients.

These are simpler than the upstream optimized rejection/table methods but
target the same distributions.

### Moments

- Univariate Fisher and Wallenius moments sum their complete finite support.
- Multivariate Fisher moments use one-dimensional marginals obtained from
  dynamic-programming normalizers.
- Multivariate Wallenius moments enumerate all feasible compositions and
  normalize the summed probabilities. As in the upstream exact-moments path,
  this can be expensive when the number of feasible compositions is very
  large.

## Validation

The included tests cover:

1. fixed univariate Fisher and Wallenius reference probabilities, CDFs,
   means, variances, and quantiles;
2. reduction of two-color multivariate distributions to the corresponding
   univariate distributions;
3. normalization of three-color Fisher and Wallenius probabilities;
4. inverse odds/urn-size approximations;
5. seeded Monte Carlo means for all four RNG families.

The source is intended to compile with standard Fortran 2018 and has no
external numerical dependencies.
