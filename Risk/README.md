# Risk for Modern Fortran

A modern Fortran 2018 and FPM translation of the computational core of the
R package `Risk` 1.0 by Saralees Nadarajah and Stephen Chan.

The library evaluates 26 financial risk measures for any continuous
probability distribution that supplies a density, CDF, quantile function,
and support bounds.

## Included algorithms

- Value at risk, expected shortfall, tail conditional median, expectiles,
  beyond value at risk, and expected proportional shortfall
- Expectation and the elementary mean-plus-standard-deviation measure
- Omega, Sortino, and Kappa measures
- Two Wang measures and two Stone measures
- Four Luce measures
- Three Sarin measures
- Four Bronshtein-Kurelenkova measures
- Scalar and vector overloads for the R routines that accept vector `alpha`
- Adaptive integration on finite, semi-infinite, and infinite intervals
- Bisection root solving

## Distribution interface

Built-in distributions are provided for:

- normal
- lognormal
- uniform
- exponential
- logistic
- Student t

The abstract `continuous_distribution` type and `callback_distribution`
allow user-defined continuous distributions. The callback example implements
a triangular distribution without changing the risk-measure code.

## Build and test

```text
fpm build
fpm test
fpm run
fpm run --example custom_distribution
```

No third-party Fortran dependency is required.

## Minimal example

```fortran
program example
   use risk
   implicit none

   type(normal_distribution) :: dist
   real(dp), parameter :: inf = huge(1.0_dp)

   dist = normal_distribution(mu=0.0_dp, sigma=1.0_dp)

   print *, varg(dist, 0.95_dp)
   print *, esg(dist, 0.95_dp)
   print *, omegag(dist, 0.0_dp, -inf, inf)
end program example
```

For vector probabilities, the same generic procedure names accept arrays:

```fortran
real(dp) :: p(3), q(3)
p = [0.90_dp, 0.95_dp, 0.99_dp]
q = varg(dist, p)
```

## Source-compatible conventions

The port follows the executable R source when its formulas differ from the
package article or from common modern terminology. In particular, `esg`
computes the lower partial mean used by the R implementation:

```text
E[X 1{X <= VaR(alpha)}] / alpha
```

This is not the usual upper-tail expected shortfall. See `PORTING.md` for all
known semantic and numerical details.

## Project layout

```text
src/       library modules
app/       demonstration program
example/   custom-distribution example
test/      automated numerical tests
original/  original R source and package metadata
LICENSES/  complete GPL license texts
```

## License

The original package declares `GPL (>= 2)`. This translation is therefore
licensed under `GPL-2.0-or-later`. Original authorship and package metadata
are retained in `NOTICE.md` and `original/`.
