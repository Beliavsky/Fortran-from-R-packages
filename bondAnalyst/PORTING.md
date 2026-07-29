# Porting notes

## Scope

The port covers all 55 computational functions exported in the original
`NAMESPACE`. The R package has no compiled source and no plotting subsystem;
its non-computational material consists mainly of roxygen/Rd documentation and
formatted presentation values. The original R source is retained under
`original/`.

## Interface changes

1. R numeric vectors become assumed-shape `real(dp)` arrays.
2. R `Date` objects become integer serial-day values. Only subtraction is
   required by the source formulas.
3. R character results produced by `format(..., scientific = FALSE)` become
   numeric `real(dp)` results after the same whole-unit rounding.
4. Functions using `polyroot` in R return the unique nonnegative economic root
   for their positive level cash-flow streams. If no such root exists, the
   Fortran function returns NaN and `ba_no_root`.
5. Every procedure has an optional `istat` argument for explicit error
   handling.

## Source behavior intentionally preserved

- `earZcbVariousPeriodicity` calculates a periodic rate, despite its name and
  effective-annual-rate terminology. The port preserves that calculation and
  adds `effectiveannualratezcb` for the conventional annual result.
- `estimatedPercentChangePVFullPrice` calculates
  `modified_duration * yield_change` without a leading minus sign. The port
  preserves it and adds `conventionalpercentpricechange`.
- `bondPriceExcessCoupon` adds the discounted excess coupon to a hard-coded
  base value of 100.
- `computingParRate` ignores the R arguments `mv` and `pv` and returns the par
  rate multiplied by 100.
- `frPricing` treats `fri` as cumulative accumulation factors and divides cash
  flows by those factors.
- `extraCompensationForHigherRisk` returns basis points by multiplying the APR
  difference by 10,000.
- The source-specific decimal rounding of each routine is retained.

## Safety corrections

The Fortran implementation rejects invalid sizes, zero day-count denominators,
non-finite values, invalid discount bases, and invalid maturity indices instead
of relying on R recycling or propagating accidental infinities.

The R duration loops use `1:(n-1)`, which behaves unexpectedly when `n=1`.
The Fortran loops correctly execute zero intermediate-coupon iterations in that
case.

The root routines use a monotone bracketed solver rather than constructing and
filtering all polynomial roots. For the positive fixed cash flows accepted by
these APIs, the economically relevant nonnegative root is unique.

## Numerical design

- No external solver or linear-algebra library is required.
- Yield and spread solvers use expanding nonnegative brackets followed by
  bisection to a relative tolerance of approximately `1e-12` before applying
  the original function's output rounding.
- Quiet NaNs are returned on errors, with status constants exposed by the
  top-level `bondanalyst` module.
