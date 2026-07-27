# Porting notes

## Type mapping

The R function's scalar arguments map to `bcc_parameters`. Its two-element
list result maps to `bcc_result`, which also carries transform probabilities,
error estimates, work counts, and status information.

The compatibility function `bcc` keeps the original argument names. The typed
`bcc_price` interface is preferred when pricing repeatedly.

## Improper integration

R's `integrate(..., upper=Inf, subdivisions=10000)` was replaced with a
self-contained adaptive Gauss-Kronrod 7/15 rule. The semi-infinite interval is
processed in successive finite panels. Convergence requires multiple recent
tail panels to be collectively smaller than the combined absolute/relative
tolerance.

This design avoids an external quadrature library and does not evaluate the
integrand exactly at `phi = 0`, where the displayed formula contains a
removable `1 / phi` singularity.

## Stable limiting calculations

The upstream examples set some long-run levels and jump intensities exactly to
zero while using tiny positive volatilities. Direct evaluation can otherwise
form `0 * NaN`. The port skips those zero-weight terms before evaluating their
intermediate logarithms. It also evaluates `1 - exp(-z)` by a short series near
zero to reduce cancellation.

These are numerical safeguards, not changes to the BCC formula.

## Complex branches

Fortran's principal complex square-root and logarithm branches are used, as in
R's standard complex arithmetic. The characteristic functions are exposed so
branch behavior can be tested independently from quadrature.

## Discounting convention

The port preserves the upstream implementation exactly: the final strike term
uses `exp(-R0 * t)`, where `R0` is the current short rate. No alternative
stochastic-bond discounting convention was substituted.

## Error handling

Invalid parameters return `status = 1` with a message. Failure to satisfy the
configured transform-tail criterion returns computed values with `status = 2`
and `converged = .false.`. Zero maturity returns intrinsic values directly.
