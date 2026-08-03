# Porting notes

## Numerical model

The backward solver follows the upstream transformed logarithmic-grid scheme.
It uses the same interior and Neumann-boundary tridiagonal coefficients,
separates survival and recovery values, and applies instrument optionality
after each implicit step. Coupon, call, put, maturity, and conversion times are
inserted into the time grid.

## Instrument representation

R reference classes were replaced by `instrument_spec`. This removes mutable
cross-instrument caches and makes repeated pricing deterministic and thread
friendly. Convertible conversion, coupon suppression after conversion, and
callable interactions cannot be represented through hidden object state; the
Fortran implementation instead applies all terms directly at each node and
time.

## Dividends

The upstream package uses R's cubic spline interpolation. The Fortran port uses
a natural cubic spline inside the grid and linear extrapolation outside it.
This is exact for linear value profiles and numerically close for smooth option
value grids, but is not bit-for-bit identical to R's default FMM spline.

## Calibration

Black-Scholes and American implied volatility use safeguarded bracketed
iterations. Intensity-link calibration is a deterministic bounded coordinate
search. It is deliberately self-contained and does not attempt to duplicate
all choices made by `limSolve`, `MASS`, or other suggested R packages.

## Term structures

A spot-rate curve is converted to piecewise-constant forward rates. A quoted
volatility curve is interpreted exactly as upstream: each quote is a total
volatility to that maturity, from which piecewise-constant forward variance is
recovered. The final forward rate and forward volatility extend indefinitely.

## Known scope boundary

The original package can price mutually dependent instrument layers using
mutable reference-class state. This port directly supports the principal
single-security cases, including a convertible bond with its own coupons,
calls, puts, recovery, and conversion. Arbitrary user-defined chains of R
objects are not emulated.

## Constructor portability

Bond constructors initialize their result records directly.  They intentionally do not forward absent optional dummy arguments through nested constructor calls, because that pattern produced optimization-sensitive behavior with GNU Fortran 14 when `-O3 -march=native` enabled aggressive inlining.  The regression suite includes a fixed convertible-bond reference value in both checked and optimized builds.
