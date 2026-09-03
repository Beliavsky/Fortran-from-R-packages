# API coverage

This document maps the supplied Rmpfr 1.1-2 computational functionality to the
native Fortran API.  "Covered" means the numerical operation is available;
it does not imply R S4 object/API compatibility.

| Upstream area | Fortran coverage | Notes |
|---|---|---|
| MPFR scalar storage / precision | Covered | `mpfr_real`, per-value precision, deep scalar assignment, finalization |
| Rounding modes | Covered | nearest, toward zero, up, down, away |
| Default precision | Covered | get/set wrappers around MPFR |
| Exponent range | Covered | get/set current range and query library limits |
| Numeric/string construction | Covered | real(dp), integer, decimal/radix string, NaN, Inf, signed zero |
| Conversion out | Covered | real(dp), integer, precision-preserving base-10 string |
| Arithmetic / comparisons | Covered | `+ - * / **`, unary minus, six comparisons |
| Elementary math | Covered | abs, sqrt, floor, ceiling, trunc, exp/expm1, logs, trig/inverse trig, hyperbolic/inverse hyperbolic |
| Gamma/special functions | Covered | gamma, lgamma, digamma, erf/erfc, zeta, Ei, Li2, Bessel J/Y 0/1/n, Airy Ai, atan2, hypot, incomplete gamma |
| Constants | Covered | pi, Euler, Catalan, log(2) |
| `sinpi`, `cospi`, `tanpi` | Covered | `tanpi` half-integer behavior deliberately corrected; see NOTICE |
| `frexpMpfr`, `ldexpMpfr` | Covered | binary decomposition/recomposition |
| `factorialMpfr` | Covered | MPFR factorial kernel |
| `chooseMpfr` | Covered | integer lower argument, including negative-n zero behavior |
| `pochMpfr` | Covered | rising factorial for nonnegative integer order |
| beta / lbeta | Mostly covered | positive arguments fully covered; upstream negative-integer special handling is retained where representable by the native integer API |
| Bernoulli numbers | Covered | upstream `-k*zeta(1-k)` convention, including B1=+1/2 |
| Summary operations | Covered | sum, product, min, max, range, cumulative sum/product; optional NaN removal for summaries |
| Matrix products | Covered | matmul, crossprod, tcrossprod in arbitrary precision |
| Normal distribution | Covered | `pnorm`, `dnorm`, and qnorm-by-inversion |
| Student-t density | Covered | central t density; upstream noncentral case was not implemented either |
| Poisson/binomial/negative-binomial densities | Covered | arbitrary-precision log/non-log kernels |
| Gamma/chi-square density, gamma CDF | Covered | incomplete-gamma backend uses MPFR |
| Integer-shape beta CDF (`pbetaI`) | Covered | finite-sum MPFR path for positive integer shapes; rational `bigq` return path omitted |
| `log1mexp`, `log1pexp` | Covered | stable branch formulas retained |
| `integrateR` | Covered | Bauer/Romberg algorithm with MPFR callbacks |
| `optimizeR` Brent | Covered | MPFR Brent minimization/maximization |
| `optimizeR` GoldenRatio | Covered | MPFR golden-section search |
| `unirootR` | Covered | Brent/zeroin-style bracketed root solver; extension policy simplified |
| `qnormI` | Covered | MPFR CDF inversion with expanding arbitrary-precision brackets |
| `hjkMpfr` | Covered | Hooke-Jeeves; deterministic coordinate order replaces randomized order |
| Rmpfr S4 classes/method dispatch | Omitted | R-specific interface layer |
| R vector recycling / `c`, `rep`, `seq`, `diff`, `apply` wrappers | Omitted | container/interface semantics rather than scalar numerical kernels |
| `mpfrArray`/S4 dimnames semantics | Omitted | native Fortran arrays are used instead |
| R `gmp::bigz` / `bigq` conversion | Omitted | R package object interoperability; GMP is still used internally by MPFR |
| R import/export serialization | Omitted | R object wire format |
| formatting/printing (`format*`, `sprintfMpfr`) | Partial | precision-preserving base-10 string export is covered; general R formatting/radix-output rules are omitted |
| plotting, vignettes, interactive helpers | Omitted | explicitly out of translation scope |

## Public modules

- `rmpfr`: umbrella module intended for users.
- `rmpfr_types`: owning MPFR scalar, operators, constructors, conversions,
  elementary/special functions, constants, precision/exponent controls.
- `rmpfr_combinatorics`: choose, Pochhammer, beta/lbeta, Bernoulli.
- `rmpfr_probability`: distribution and stable-log kernels.
- `rmpfr_summary`: scalar summaries and cumulative output routines.
- `rmpfr_matrix`: arbitrary-precision matrix products.
- `rmpfr_algorithms`: integration, root finding, qnorm, scalar optimization,
  and Hooke-Jeeves.
