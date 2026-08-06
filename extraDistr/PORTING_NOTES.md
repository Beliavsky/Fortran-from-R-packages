# Porting notes

## Translation strategy

The upstream package implements most numerical kernels in C++ through Rcpp.
The Fortran port replaces Rcpp/Rmath calls with self-contained implementations
of the normal, beta, gamma, Student-t, binomial, negative-binomial, Poisson,
and Bessel functions and their inverses where needed.

The public `d/p/q/r` names are preserved. Fortran scalar procedures replace R
vector recycling. Random procedures accept a scalar number of draws and return
allocatable vectors or matrices.

## Random-number generation

`seed_rng` deterministically initializes the compiler's intrinsic random-number
state. Normal draws use Box-Muller; gamma draws use Marsaglia-Tsang; other
families use transformations, mixtures, or inverse-CDF generation. A fixed seed
is reproducible within the same Fortran implementation, but streams are not
bit-for-bit identical to R's RNG.

## Numerical methods

- Normal quantiles use a rational approximation with Newton refinement.
- Incomplete beta and gamma functions use continued fractions and series.
- General quantiles use safeguarded bisection where no closed form exists.
- Skellam probabilities use an integer-order modified-Bessel recurrence.
- Mixture probabilities are normalized internally from nonnegative weights.

## Source-compatible edge behavior

The upstream C++ has three observable edge cases that are not the conventional
mathematical definitions. They remain the default:

1. At `x == mu`, `dslash` omits division by `sigma`.
2. `pdlaplace` chooses its CDF branch from `q < 0`, not `q-location < 0`.
3. `ddirichlet` does not test whether components sum to one.

Each procedure has an optional `source_compatible` argument. Set it to
`.false.` for the conventional behavior.

## Deliberate interface differences

- Discrete variates are represented by integers in density/CDF/quantile APIs,
  so the R behavior of returning zero for a noninteger input is represented by
  the Fortran type system rather than a runtime branch.
- `rdlaplace` returns real values because upstream permits a noninteger
  location shift.
- `ddgamma`/`pdgamma`/`rdgamma` and `dgpois`/`pgpois`/`rgpois` preserve
  the R-style `rate` argument and also accept an optional overriding `scale`.
- Categorical labels are not part of the numerical library; returned values are
  one-based integer category indices.
- Invalid scalar parameters generally produce IEEE NaN, analogous to upstream
  warnings plus NaN. Shape mismatches are represented by zero/NaN results as
  documented by each procedure rather than R exceptions.
