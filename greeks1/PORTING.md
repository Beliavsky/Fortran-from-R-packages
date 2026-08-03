# Porting notes

## R interfaces

R lists, named vectors, tibbles, S3 dispatch, and `...` are represented with
explicit Fortran arguments and the `greek_result` derived type. Input
vectorization is intentionally explicit: callers loop over parameter vectors.

The R function `Greeks` is named `option_greeks` in Fortran to avoid a name
collision with the public module `greeks`.

## Random numbers

The R package uses `dqrng`. This translation uses a deterministic Park-Miller
uniform generator and Box-Muller normal generation, so identical seeds do not
produce identical R and Fortran paths. The statistical model and estimators are
unchanged. Antithetic generation is supported for all translated Monte Carlo
routines.

The upstream jump-diffusion default draws Student t(3) jump sizes. The Fortran
translation includes that default. An arbitrary R jump-generator closure is
not directly portable; changing the jump law requires adapting the simulation
routine.

## Custom payoffs

The general Malliavin routines accept optional pure Fortran callbacks. The
arithmetic-Asian derivative estimators also accept a derivative callback. The
control-variate estimator remains restricted to calls and puts because it
requires an exact geometric-Asian control with a known closed form.

## American tree

The original R wrapper computes a corrected price:

`American tree + exact European - European tree`.

The Fortran routine preserves this behavior and then computes Greeks by finite
differences. The raw tree values are separately available through
`binomial_values`.

## Implied volatility

The exact European routine preserves the upstream Halley update. The general
routine uses deterministic repeated evaluations and Newton updates. For Monte
Carlo option types, convergence depends on path count, seed, and tolerance.

## Numerical differences and fixes

- The exact geometric-Asian put vomma is ported from the explicit upstream C++
  branch; it is not inferred by call/put symmetry.
- Array bounds and dimensions are checked before simulation.
- The path constructor avoids the upstream inclusive-range ambiguity in
  `make_BM` by consuming exactly `paths * steps` increments.
- Monte Carlo results include estimated standard errors.
- Status codes replace R exceptions in the numerical library.
