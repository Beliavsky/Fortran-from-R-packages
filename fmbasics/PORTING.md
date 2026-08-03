# Porting notes

## Object model

The R package uses S3 classes and method dispatch. The Fortran port uses
strongly typed derived types and explicit procedures. For example, R's generic
`interpolate()` maps to `interpolate_zero`, `interpolate_credit` or
`interpolate_vol`, depending on the object being evaluated.

R vector recycling is implemented where it is natural for rate and discount
objects. Fortran callers otherwise pass conformable arrays, making shape
requirements explicit.

## Dates and calendars

The upstream package imports date calculations from `fmdates`. This project
stores a civil date as an integer day count relative to 1970-01-01 and includes
its own Gregorian conversion, day-count and business-day routines.

Named calendars provide weekends and selected recurring/common holidays needed
by the package's benchmark conventions and tests. They are deterministic and
portable, but are not a full historical or future market-holiday database.
One-off holidays, local/state holidays and rule changes may differ from
`fmdates`. Production settlement systems should supply or validate their own
calendar data.

## Credit bootstrap

The original `CDSCurve` conversion delegates to `credule`. The Fortran port
uses a self-contained bootstrap with piecewise-constant hazard rates, quarterly
premium periods by default, accrual-on-default approximation, supplied recovery
rate and discount factors from `zero_curve_t`. This preserves the computational
purpose without an external dependency, but exact results can differ from
`credule` because of schedule and accrual conventions.

## Interpolation

- Constant, linear and natural cubic one-dimensional interpolation are
  implemented internally.
- Log-discount-factor zero-curve interpolation is linear in log discount
  factors and flat in the terminal zero rates outside the pillar range.
- Volatility interpolation is linear in total variance through maturity and
  natural cubic through smile coordinates.

## Operators and presentation

R's overloaded arithmetic, comparisons, concatenation and subsetting are
represented by named functions, constructors and direct component access.
Formatting, printing, `as_tibble`, `tbl_sum`, `type_sum`, validation predicates
whose only role is R class checking, and package-load hooks are omitted.

## Input data

The original example zero-curve and volatility-surface CSV files are retained
under `data/`. CSV readers are intentionally small and expect the same layout as
those files; they are not general quoted-field CSV parsers.
