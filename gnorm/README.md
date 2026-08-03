# gnorm-fortran

Modern Fortran translation of the computational code in the R package
`gnorm` 1.0.2.

The library implements the generalized normal distribution, also called the
exponential-power distribution, with density

```text
f(x) = beta * exp(-(|x-mu|/alpha)^beta) /
       (2 * alpha * Gamma(1/beta)).
```

## Implemented API

- `dgnorm` - density or log density
- `pgnorm` - lower- or upper-tail distribution function
- `qgnorm` - lower- or upper-tail quantile
- `rgnorm` - allocate and return random deviates
- `rgnorm_fill` - fill a caller-provided array with random deviates
- `gnorm_mean` and `gnorm_variance` - theoretical moments

`dgnorm`, `pgnorm`, and `qgnorm` are elemental functions, so scalar calls work
unchanged with conformable arrays.

## Build with FPM

```text
fpm build
fpm test
fpm run --example density_cdf_example
fpm run --example quantile_example
fpm run demo_gnorm
```

The project has no external dependencies.

## Minimal example

```fortran
use gnorm

real(dp) :: probability, quantile
real(dp), allocatable :: draws(:)

probability = pgnorm(1.0_dp, alpha=sqrt(2.0_dp), beta=2.0_dp)
quantile = qgnorm(0.975_dp, alpha=sqrt(2.0_dp), beta=2.0_dp)
draws = rgnorm(1000, mu=0.0_dp, alpha=1.0_dp, beta=1.5_dp, seed=42_i8)
```

The choice `alpha=sqrt(2)` and `beta=2` gives the standard normal
distribution. `beta=1` gives a Laplace distribution with scale `alpha`.

## Porting notes

The regularized incomplete gamma function, its inverse, and gamma random
sampling are implemented in standard Fortran. Invalid nonpositive `alpha` or
`beta` values return IEEE NaNs for distribution functions and an explicit
status for random generation.

The upstream `qgnorm` branch for `lower.tail=FALSE, log.p=TRUE` applies
`log(1-p)` to a log probability. The Fortran implementation uses the intended
R distribution-function semantics, namely `1-exp(p)`.

See `API.md`, `PORTING.md`, and `TRANSLATION_COVERAGE.md` for details.

## License

GPL-2.0-or-later. Original package sources and metadata are retained under
`original/`.
