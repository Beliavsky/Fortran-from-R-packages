# Porting notes

## Scope

All package-specific numerical work was translated. Graphics, S3 dispatch,
R data frames, list/closure construction, and helper string operators were
not translated because they are R runtime infrastructure rather than yield
curve algorithms.

## Wilson function

The R source evaluates

```text
exp(-ufr*(t+u)) *
(alpha*min(t,u) - sinh(alpha*min(t,u))/exp(alpha*max(t,u)))
```

The Fortran code evaluates the algebraically identical hyperbolic term as

```text
0.5 * (exp(-alpha*(max-min)) - exp(-alpha*(max+min)))
```

This avoids intermediate overflow when `alpha*min(t,u)` is large.

## Calibration solver

R uses `solve(C W C^T)`. The translation uses Gaussian elimination with
partial pivoting. A singular or numerically rank-deficient calibration matrix
returns `sw_singular_system`; no implicit ridge regularization is added.

## Instrument conventions

The upstream simple conventions are preserved:

- A LIBOR instrument has one payment at tenor with cashflow
  `1 + rate*tenor`.
- A swap pays `rate/frequency` at each date and returns notional at maturity.
- A bond uses the same coupon convention and returns notional at maturity.
- A bond schedule is generated backward from maturity, so an off-cycle tenor
  produces a short first interval but still receives a full regular coupon,
  exactly as in the R code.
- LIBOR and swap market values are one; bond market value is its price.

The package intentionally does not implement business-day calendars, day
count conventions, floating-leg projections, or accrued interest.

## Input validation

The upstream `sequence(tenor*frequency)` operation only works sensibly for a
positive integer payment count. The Fortran API checks that condition
explicitly and reports a diagnostic rather than relying on an R coercion or
runtime error.

## Curve representation

R returns closures for the pricing and compound-kernel functions. Fortran
returns `type(smith_wilson_curve)`, which owns the calibration data and
provides type-bound scalar and vector methods. This avoids procedure-lifetime
issues and makes copying and serialization straightforward.

## Precision

All public real quantities use

```fortran
integer, parameter :: dp = kind(1.0d0)
```

No `kind=8` assumptions are used.
