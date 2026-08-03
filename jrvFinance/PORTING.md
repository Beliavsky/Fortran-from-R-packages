# Porting notes

## R-to-Fortran mappings

R function names containing dots use underscores in Fortran. R lists and Date
objects are replaced by strongly typed derived types. Vectorized bond wrappers
accept arrays of equal length; unlike R `mapply`, they do not silently recycle
shorter arrays.

Continuous compounding is encoded as frequency `0.0_dp`, avoiding a dependency
on constructing an IEEE infinity at each call.

## Numerical behavior

The formulas for NPV, IRR, annuities, coupon schedules, accrued interest,
Black-Scholes prices, and Greeks follow the upstream R source. The root layer
uses safeguarded Newton steps and a deterministic geometric bracketing fallback.
Implied volatility uses the upstream Brenner-Subrahmanyam-style initial guess.

## Deliberate corrections

The upstream short-bond duration branch calls `yearFraction` without its
required reference dates. The Fortran implementation returns actual days divided
by 365 for a single-payment short instrument.

Invalid inputs and failed numerical solves return explicit status codes and,
where a scalar result is required, an IEEE quiet NaN. They do not print R-style
warnings or terminate the process.

## Omitted infrastructure

R documentation rendering, vignettes, dynamic argument recycling, and R object
printing are not part of the compiled library. Original R sources and manuals
are retained under `original/`.
