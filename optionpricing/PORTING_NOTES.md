# Porting notes

## Data model

R named vectors and matrices are represented by typed results:

- `european_result`: price, delta, correct gamma, and upstream gamma
- `greeks_result`: three estimates, three 95 percent half-widths, and status
- `moments_result`: deterministic price, delta, and gamma quantities
- `conditional_result`: conditional means and six control-variate means

Arrays follow normal Fortran column-major conventions.

## Black-Scholes gamma correction

Upstream uses the normal CDF in the gamma numerator. This is not the
Black-Scholes gamma formula. The port uses the normal density and places the
upstream calculation in `upstream_gamma` so that the correction is explicit
and testable.

## Regression

R's `lm` is replaced by rank-revealing, column-pivoted modified Gram-Schmidt QR
with reorthogonalization. This matters because the six control variables can
be nearly dependent in small pilot samples. Rank-deficient columns receive
zero coefficients rather than destabilizing the estimate.

## Lord integration

R integrates from negative infinity. The adaptive Fortran integrator starts at
the smaller of -12 and eight normal standard deviations below the upper bound.
The omitted normal-tail contribution is far below the integration tolerance
for ordinary double precision.

## Random numbers

The stochastic estimators accept optional integer seeds. Repeated calls with
the same compiler and seed are reproducible. Fortran does not standardize the
exact `random_number` sequence across compilers, so statistical tests rely on
invariants and error bounds rather than compiler-specific random streams.

## Randomized QMC

The functional QMC method in the R package is a randomly shifted Korobov
lattice, optionally followed by the Baker transform. It is fully translated.
The source's Sobol branch only prints a message and comments out the external
Sobol call, so the port does not advertise Sobol as an available method.

The transformation modes are:

- `std`: conditional standard-normal projection
- `pca`: principal-component factorization
- `pcamain`: leading PCA directions completed to an orthonormal basis
- `lt`: linear transformation directions from payoff sensitivities
- `ltpca`: LT rotation applied after conditional PCA

## Error estimates

Monte Carlo half-widths use 1.96 times the sample standard deviation divided
by the square root of the number of production observations. QMC half-widths
use variation across independently randomized lattice copies, matching the R
package's design.

## Source conventions

- Fortran 2018 free form
- Lower-case source
- `implicit none` in every module and program
- `dp = kind(1.0d0)`
- No line longer than 132 characters
- ASCII-only release text

## Cross-compiler deterministic references

The fixed conditional-Monte-Carlo regression test evaluates exponentials,
`erfc`, cancellation-prone probability differences, and a Newton solve. The
last few floating-point bits therefore depend on the compiler and math library.
The test uses combined absolute and relative tolerances that still require
about 10 to 11 significant decimal digits. The numerical implementation and
published reference values are unchanged.
