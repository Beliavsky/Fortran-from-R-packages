# Porting notes

## Source and scope

The source archive contains `tvm` 0.5.2, an R package for time-value-of-money,
loan, and rate-curve calculations. It declares `MIT + file LICENSE` and names
Juan Manuel Truppia as author and copyright holder.

All numerical routines were translated. Plotting and R-specific object,
documentation, and testing infrastructure were not translated into compiled
Fortran.

## Data representation

R numeric vectors become `real(dp)` arrays, where
`dp = kind(1.0d0)`. Loans and curves become the typed objects `loan_t` and
`rate_curve_t` rather than untyped R lists with S3 classes.

The original `xnpv` and `xirr` accept R `Date` values or year fractions. The
Fortran generic interfaces accept either integer day serials or real year
fractions. Day serials need not use any particular calendar epoch because only
differences from the first date are used.

## Root finding

R's `uniroot` calls are replaced by a self-contained bracketed bisection
solver. IRR routines extend the upper end of the search interval when needed,
mirroring the intent of `extendInt = "yes"`. Optional status outputs distinguish
successful and unsuccessful searches.

Financial cashflows can have multiple IRRs. As in the R package, the returned
root depends on the supplied search interval.

## Curve interpolation

The R constructor defaults to `splinefun(..., method = "monoH.FC")`. The
Fortran implementation uses a Fritsch-Carlson monotone piecewise-cubic Hermite
interpolator, including monotonicity-preserving endpoint slopes. It preserves
shape and node values but is not claimed to be bit-for-bit identical to every
R version's spline implementation.

A direct discount-factor callback is retained and evaluated directly between
knots. A rate callback is sampled at the supplied knots before bootstrapping,
matching the R constructor's structure.

## Rate scaling

The original classification is preserved:

- nominal scaling: `zero_nom`, `german`, `french`, `swap`, `fut`, `zero_cont`;
- effective scaling: `zero_eff`.

Although continuously compounded rates are not usually described as nominal,
this classification intentionally follows the upstream `unscale` and `rescale`
functions.

## Safeguards and corrections

The translation adds the following defined behavior and validation:

- `pmt` returns `amount / maturity` at a zero rate instead of evaluating `0/0`;
- invalid dimensions, non-increasing knots, unsupported rate types, and invalid
  grace periods are rejected explicitly;
- discount transformations that require positive factors validate them;
- requested grace-period metadata is retained in `loan_t`; the R constructor
  stores zero in those metadata fields even when nonzero grace arguments are
  supplied, although its cashflow calculation uses the arguments;
- one-knot output rate curves are handled as constant-rate curves.

These changes do not alter the ordinary positive-rate examples and tests in the
original package.

## Procedure callbacks

Fortran callbacks must match the scalar interface

```fortran
function f(x) result(y)
   use tvm_kinds, only : dp
   real(dp), intent(in) :: x
   real(dp) :: y
end function f
```

Callbacks can be module or external procedures. A curve that stores a discount
callback should not outlive the procedure it references.

## Presentation

`plot.rate_curve` is intentionally omitted. `rate_curve_t%rates`,
`rate_curve_t%rate_grid`, and `rate_curve_t%discount` return the complete data
needed for plotting or export.
