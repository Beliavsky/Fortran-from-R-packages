# nbconv-fortran

Modern Fortran translation of the computational code in the R package
`nbconv` 1.0.1 by Gregory Bedwell.

The package evaluates convolutions of independent negative-binomial random
variables with different means/probabilities and dispersion parameters.
It implements all computational methods exported by the R package:

- Furman's series representation (`exact`)
- moment-matched negative-binomial approximation (`moments`)
- saddlepoint approximation (`saddlepoint`)
- PMF, CDF, quantile, and random-generation interfaces
- mean, variance, skewness, excess kurtosis, and Furman K-distribution mean

The translation is self-contained and has no external numerical-library
dependencies.

## Build

With FPM:

```text
fpm build
fpm test
fpm run --example basic
```

Or compile the files in `src/` in this dependency order:

```text
nbconv_kinds.f90
nbconv_rng.f90
nbconv_math.f90
nbconv_exact.f90
nbconv_approximations.f90
nbconv_api.f90
nbconv.f90
```

## Basic use

```fortran
use nbconv, only : dp, dnbconv, qnbconv, nbconv_summary, nbconv_params

real(dp), parameter :: mus(2) = [100.0_dp, 10.0_dp]
real(dp), parameter :: phis(2) = [5.0_dp, 8.0_dp]
integer, parameter :: counts(3) = [50, 100, 150]
real(dp), allocatable :: pmf(:)
type(nbconv_summary) :: s

pmf = dnbconv(counts, mus, phis, parameterization="mu", method="exact")
s = nbconv_params(mus, phis, parameterization="mu")
```

For probability-of-success parameterization, pass `parameterization="p"`.
The explicit routines `dnbconv_mu`, `dnbconv_p`, `pnbconv_mu`,
`pnbconv_p`, and so on are also available and avoid the runtime
parameterization string.

## Public API

Low-level evaluators:

- `nb_sum_exact(ps, phis, counts, ...)`
- `nb_sum_moments(mus, phis, counts)`
- `nb_sum_saddlepoint(mus, phis, counts, normalize)`

R-style distribution operations:

- `dnbconv(...)`
- `pnbconv(...)`
- `qnbconv(...)`
- `rnbconv(...)`
- `nbconv_params(...)`

Type-safe parameterization-specific variants append `_mu` or `_p`.
`nbconv_seed()` seeds the package's reproducible Park-Miller random-number
stream.

## License

The upstream R package is licensed GPL (>= 3). This translation is therefore
released under GPL-3.0-or-later. See `LICENSE` and `UPSTREAM.md`.
