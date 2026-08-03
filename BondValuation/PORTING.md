# Porting notes

## Source version and license

This project translates the computational code of R package `BondValuation`
0.1.1, published 2022-05-28. The original package declares `GPL-3`; this port
uses SPDX identifier `GPL-3.0-only` and retains the full original source tree.

## Algorithm mapping

| R/native component | Fortran component |
|---|---|
| `AccrInt` | `accrued_interest`, `accr_int`, compatibility `AccrInt` |
| `AnnivDates` | `build_bond_schedule`, `anniv_dates`, compatibility `AnnivDates` |
| `DP` | `dirty_price`, `dp_value`, compatibility `DP` |
| `BondVal.Price` | `bond_price`, `bond_val_price`, compatibility `BondVal_Price` |
| `BondVal.Yield` | `bond_yield`, `bond_val_yield`, compatibility `BondVal_Yield` |
| C++ date kernels | `bondvaluation_dates` and `bondvaluation_compat` |
| `DIST`/`PayCalc` | `day_count_fraction`/`PayCalc` |
| price derivatives | `dm_MyPriceEqn`, `ModDUR`, `CONV` |
| `NonBusDays.Brazil` | exact embedded bitset and text provenance table |

## Preserved source behavior

- Coupon rates and yields are supplied in percent.
- Coupon frequencies are restricted to the package's supported values.
- Odd first and last coupons are valued from day-count fractions.
- End-of-month behavior is inferred from supplied coupon anchors by default.
- Business/252 uses the package's exact Brazilian non-business-day table.
- The source's Actual(No Leap)/365 behavior subtracts at most one leap day from
  an interval, and is intentionally preserved.
- Convexity includes the source's factor of one half. It is therefore the
  coefficient convention used by the package rather than the unscaled second
  derivative sometimes reported by fixed-income libraries.
- A final remaining payment may be discounted by the source's simple-period
  convention; this is controlled by `simple_last_period`.

## Intentional interface changes

- R `Date`, `timeDate`, data frames, vectors with recycling, and `NA` are replaced
  by typed Fortran dates, terms, schedules, results, warning flags, and status codes.
- The scalar numerical API is explicit. Large bond panels are processed with a
  normal Fortran loop, avoiding R-style implicit recycling and coercion.
- R warnings that merely describe inferred/dropped dates are represented by
  `schedule_warnings` fields.
- Yield solving uses safeguarded bisection instead of relying solely on the
  original Newton iteration. The translated Newton and analytical derivative
  kernels remain available in `bondvaluation_compat`.
- The original approximate previous/successor date kernels used rounded day
  offsets before repairing day-of-month. The preferred schedule engine uses
  exact calendar-month arithmetic; source-compatible kernels are still exposed.

## Corrections and robustness improvements

- Invalid dates, orderings, coupon frequencies, settlements, and discount bases
  return deterministic status codes rather than relying on R coercion failures.
- Coupon-date candidates are sorted and deduplicated before use.
- Yield roots are bracketed adaptively and cannot silently diverge.
- The last-payment direct yield formula stores its value before assigning the
  complete result, avoiding alias-sensitive component assignment.
- Invalid Newton probes are rejected using IEEE finite-value checks.

## Not translated as computational APIs

- Rcpp registration and R namespace machinery
- R data-frame construction and display formatting
- `timeDate` class metadata
- bundled demonstration tables as R objects

The original datasets and documentation remain under `original_source/`.
