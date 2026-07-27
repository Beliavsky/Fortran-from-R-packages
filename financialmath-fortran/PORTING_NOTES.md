# Porting notes

## R `NA` arguments

The original package uses `NA` to indicate which scalar should be solved for.
Fortran uses explicit strings such as `unknown='rate'` or
`unknown='payment'`. This removes ambiguity and avoids using NaN as ordinary
control flow.

## Root finding

Upstream annuity and IRR routines often construct a polynomial and call R's
`polyroot`, then filter rounded roots. The Fortran port evaluates the financial
formula directly and uses bounded bisection. This avoids complex intermediate
roots and is more portable, but it returns one real root from the requested
bracket rather than every polynomial root.

`irr` defaults to a lower bound near -100 percent and expands its positive upper
bound when needed. Optional bounds allow callers to isolate another root when a
non-conventional cash-flow stream has multiple IRRs.

## Rate conventions

- `solve_tvm` interprets `rate` as a nominal annual interest rate convertible
  `compounding_frequency` times per year.
- Annuity and loan routines convert that nominal rate to the payment frequency.
- Forward routines retain the upstream continuously compounded `r` convention.
- `bls_order1` retains the upstream continuous dividend-yield convention.

## Amortization

For an unknown term, the upstream package reports a fractional last period and
uses an approximate partial-period formula. The Fortran schedule accrues the
remaining balance by `(1+j)^fraction-1` and then pays it off exactly. Balloon
and drop-payment diagnostics are still returned.

## Payoff tables

The original plotting functions construct rounded grids mainly for display.
The Fortran routines return deterministic 11-point grids and preserve the
financial payoff/profit formulas. Callers that require a custom grid can apply
the same scalar payoff formulas directly to their own terminal-price array.

## Numerical robustness

- Zero-maturity and zero-volatility Black-Scholes prices use deterministic
  limits in the scalar call/put routines.
- Root-solving routines return a typed status instead of terminating the
  process.
- Opposite option positions and bull/bear spreads are tested for exact payoff
  and profit cancellation.
- All array dimensions are checked by the Fortran runtime in the validation
  build.

## Linker note

GNU Fortran may emit an executable-stack warning on Linux because internal
procedure callbacks are used for root solving. This is a compiler trampoline
implementation detail; it does not affect the tested numerical results or FPM
layout.
