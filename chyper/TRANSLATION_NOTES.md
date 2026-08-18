# Translation notes

## Upstream

- R package: `chyper` 0.3.1
- Author: William Nickols
- License: MIT
- Upstream computational source: `R/chyper.R`

The unmodified upstream package snapshot is retained under `upstream/` for
provenance.

## Numerical translation

The R implementation evaluates the overlap law recursively, using repeated
calls to R's `dhyper()` and a recurrence intended to avoid probabilities that
round to zero. The Fortran translation expresses the same stochastic recursion
directly:

1. For the first sample, the number selected from the common region is an
   ordinary hypergeometric variate.
2. Conditional on `x` common objects surviving all previous samples, the next
   sample keeps `y` of them according to a hypergeometric distribution with
   `x` successes and `n(i)+s-x` failures.
3. Dynamic programming propagates the probability mass through all samples.

Hypergeometric probabilities are computed from `log_gamma`, which avoids the
large intermediate binomial coefficients that would overflow integer types.
The resulting probability vector is normalized once at the end to remove tiny
floating-point accumulation drift.

This implementation has the same distributional target as the upstream nested
R recursion but is simpler, fully bounds-safe, and does not require the
round-to-zero detection logic.

## Quantile convention

The upstream `qchyper()` returns the first integer for which the CDF is
**strictly greater** than `p` (`which(cumsum > p)`). The Fortran `qchyper`
preserves that convention rather than silently switching to `CDF >= p`.

## MLE routines

`mle_s`, `mle_n`, and `mle_m` retain the upstream integer hill-search strategy
based on the stated unimodality of the likelihood. The lower bound for `mle_s`
is written as the mathematically necessary `max(0,max(m-n))`; this avoids
asking the PMF evaluator to score invalid populations.

For `mle_n`, upstream returns R `Inf` when every observed overlap is zero. The
Fortran routine has a real result and returns `huge(1.0_dp)` for that boundary
case.

## Omitted R-specific code

There is no plotting code in the upstream package. R-specific vector recycling,
messages, data frames used only for sorting, roxygen documentation, and dynamic
R error handling are not reproduced. Explicit scalar and vector Fortran entry
points are provided instead.

## v0.1.1 test portability fix

`test_quantile_pvalue` in v0.1.0 used the literal
`0.9478252815000214` and expected quantile 2.  High-precision evaluation gives
`F(1) = 0.94782528150002152109...`, so that literal is actually below the CDF
jump by about `1.21e-16`; with the upstream strict `CDF > p` convention the
correct quantile is therefore 1.  Some compilers rounded the computed CDF down
to the same double as the literal and returned 2, making the test nonportable.
The v0.1.1 test uses `p=0.95`, safely inside the `q=2` interval.  No library
algorithm or public API behavior was changed.

