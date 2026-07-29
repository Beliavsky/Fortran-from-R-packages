# bondAnalyst-fortran

A modern Fortran/FPM translation of the computational routines in the R package
**bondAnalyst 1.0.1**, by MaheshP Kumar.

The library implements all 55 functions exported by the original package:

- coupon-bond, zero-coupon bond, FRN, Treasury-bill, commercial-paper, and
  money-market valuation;
- yield-to-maturity, yield-to-call, discount-margin, G-spread, and Z-spread
  calculations;
- spot-rate, par-rate, and implied-forward-rate calculations;
- accrued interest, Macaulay duration, modified duration, effective duration,
  money duration, and PVBP;
- APR periodicity conversion and zero-coupon yield conversions.

The original R names are retained as lower-case Fortran procedure names.
Fortran is case-insensitive, so a call such as `computingZspread` is equivalent
to `computingzspread`.

## Requirements

- A Fortran 2018 compiler
- FPM for normal package builds

The library has no external numerical-library dependency.

## Build and test

```text
fpm build
fpm test
fpm run
fpm run --example curve_and_duration
```

## Minimal example

```fortran
program example
   use bondanalyst, only : dp, bondpriceyearlycoupons
   implicit none

   real(dp), parameter :: coupons(5) = [4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp]
   real(dp), parameter :: times(5) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]

   print *, bondpriceyearlycoupons(coupons, times, 100.0_dp, 5, 0.06_dp)
end program example
```

This prints `91.575`, matching the source package's documented example and
three-decimal rounding convention.

## Error handling

Every public procedure accepts an optional final integer `istat` argument.
A successful call sets it to `ba_success`. Invalid calls return a quiet NaN and
one of:

- `ba_invalid_argument`
- `ba_size_mismatch`
- `ba_no_root`
- `ba_out_of_range`

Omitting `istat` gives a compact function-call interface similar to R.

## Important compatibility notes

The port intentionally preserves source-code behavior, even where a name or
formula is unconventional. In particular:

- `earzcbvariousperiodicity` returns the periodic rate calculated by the R
  source, not an effective annual rate. The added
  `effectiveannualratezcb` function gives the conventional effective annual
  rate.
- `estimatedpercentchangepvfullprice` preserves the R source's sign. The added
  `conventionalpercentpricechange` uses `-modified_duration * yield_change`.
- `computingparrate` returns a percentage value and retains the source
  signature, although its `mv` and `pv` inputs do not enter the formula.
- `bondpriceexcesscoupon` retains the source's fixed base price of 100.

See [PORTING.md](PORTING.md) for the complete mapping and differences.

## Provenance and license

The original R source files, `DESCRIPTION`, `NAMESPACE`, and references are
retained under `original/` for review and provenance.

The original package declares `License: GPL-3`. This port is distributed under
**GPL-3.0-only** and includes the complete license text in [LICENSE](LICENSE).
