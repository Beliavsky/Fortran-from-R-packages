# Translation notes

## Upstream

- Package: `LindleyPowerSeries`
- Version: 1.0.1
- Authors: Saralees Nadarajah, Yuancheng Si, Peihao Wang
- Upstream license: GPL (>= 2)
- CRAN DESCRIPTION reference DOI: 10.1007/s13171-018-0150-x
- Upstream language: pure R; `NeedsCompilation: no`
- Imported runtime packages: `stats`, `lamW (>= 1.3.0)`

The original `DESCRIPTION`, `NAMESPACE`, and R sources are retained under
`upstream/`.

## Function coverage

| R function | Fortran function |
|---|---|
| `plindleygeometric` | `plindleygeometric` |
| `dlindleygeometric` | `dlindleygeometric` |
| `hlindleygeometric` | `hlindleygeometric` |
| `qlindleygeometric` | `qlindleygeometric` |
| `rlindleygeometric` | `rlindleygeometric` |
| `plindleylogarithmic` | `plindleylogarithmic` |
| `dlindleylogarithmic` | `dlindleylogarithmic` |
| `hlindleylogarithmic` | `hlindleylogarithmic` |
| `qlindleylogarithmic` | `qlindleylogarithmic` |
| `rlindleylogarithmic` | `rlindleylogarithmic` |
| `plindleynb` | `plindleynb` |
| `dlindleynb` | `dlindleynb` |
| `hlindleynb` | `hlindleynb` |
| `qlindleynb` | `qlindleynb` |
| `rlindleynb` | `rlindleynb` |
| `plindleybinomial` | `plindleybinomial` |
| `dlindleybinomial` | `dlindleybinomial` |
| `hlindleybinomial` | `hlindleybinomial` |
| `qlindleybinomial` | `qlindleybinomial` |
| `rlindleybinomial` | `rlindleybinomial` |
| `plindleypoisson` | `plindleypoisson` |
| `dlindleypoisson` | `dlindleypoisson` |
| `hlindleypoisson` | `hlindleypoisson` |
| `qlindleypoisson` | `qlindleypoisson` |
| `rlindleypoisson` | `rlindleypoisson` |

There is no plotting code in the upstream package.

## Mathematical structure

Let

`G(x) = 1 - ((lambda + 1 + lambda*x)/(lambda + 1))*exp(-lambda*x)`

be the ordinary Lindley CDF and let `g(x)` be its density. For a power-series
normalizer `A`, all five families satisfy

`F(x) = A(theta*G(x))/A(theta)`

and

`f(x) = theta*g(x)*A'(theta*G(x))/A(theta)`.

The five choices are:

- geometric: `A(t) = t/(1-t)`;
- logarithmic: `A(t) = -log(1-t)`;
- negative binomial: `A(t) = (t/(1-t))**m`;
- binomial: `A(t) = (1+t)**m - 1`;
- Poisson: `A(t) = exp(t) - 1`.

This common representation was used to simplify and cross-check the R
expressions.

## Corrected upstream defects

Two source expressions are inconsistent with the package's own documented
power-series definition and fail basic distribution identities. They are
corrected in the Fortran translation:

1. `qlindleybinomial`: upstream computes
   `(p*A(theta) - 1)^(1/m) - 1`. Since
   `A(t)=(1+t)^m-1`, the inverse is
   `A^{-1}(y)=(1+y)^(1/m)-1`; the sign must therefore be `+1`.
   The upstream quantile/hazard/RNG wrappers also restrict `theta < 1`, while
   the documentation and CDF/PDF correctly specify `theta > 0`. The Fortran
   implementation supports all positive `theta`.

2. `hlindleypoisson`: upstream multiplies the correct hazard expression by an
   additional `theta*G(x)` factor. The Fortran implementation uses the general
   identity `h(x)=f(x)/(1-F(x))`, algebraically evaluated as
   `theta*g(x)/(exp(theta*(1-G(x)))-1)`.

These fixes are covered by explicit tests.

## Interface differences

- R vector arguments become scalar `elemental` Fortran functions, which may be
  called on arrays directly.
- The RNG entry points take an integer `n` and return an allocatable rank-1
  array.
- The negative-binomial and binomial `m` argument is an integer in Fortran,
  matching the intended positive-integer use.
- Invalid-argument checking is not reproduced as R `stopifnot` calls inside
  every elemental function. Callers should obey the documented parameter
  domains.

## Validation

The test suite includes:

- high-precision fixed PDF/CDF/hazard reference values for all five families;
- high-precision quantile references and CDF/quantile round trips;
- the hazard identity `h=f/(1-F)`;
- `log_p` behavior;
- `theta > 1` Lindley-binomial quantiles;
- large-parameter Poisson/binomial stability;
- far-tail hazard stability;
- RNG size/support checks.

The release was compiled with GNU Fortran 14.2 using
`-std=f2018 -Wall -Wextra -Wimplicit-interface -fcheck=all`.
