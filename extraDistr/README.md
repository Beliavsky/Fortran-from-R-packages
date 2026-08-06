# extraDistr-fortran

Modern Fortran translation of the computational routines in R package
`extraDistr` 1.10.0.5.

The library provides the complete 194-name numerical API exported by the R
package: densities or probability masses, cumulative probabilities, quantiles,
random generation, mixtures, and multivariate distributions. It is
self-contained and requires no BLAS, LAPACK, R, Rcpp, or RcppArmadillo runtime.

## Main modules

- `extra_distr`: primary public module
- `extra_distr_continuous`: continuous univariate distributions
- `extra_distr_discrete`: discrete univariate distributions
- `extra_distr_multivariate`: mixtures and multivariate distributions
- `extra_distr_math`: probability and special-function kernel
- `extra_distr_rng`: seeded random-number utilities

All floating-point calculations use

```fortran
integer, parameter :: dp = kind(1.0d0)
```

## Example

```fortran
program example
  use extra_distr
  implicit none
  real(dp) :: x
  real(dp), allocatable :: draws(:)

  call seed_rng(12345)
  x = qgpd(0.99_dp, mu=0.0_dp, sigma=1.5_dp, xi=0.2_dp)
  draws = rlaplace(1000, mu=0.0_dp, sigma=1.0_dp)

  write(*,'(a,f12.6)') 'GPD quantile: ', x
  write(*,'(a,f12.6)') 'Laplace mean: ', sum(draws)/real(size(draws),dp)
end program example
```

Scalar `d`, `p`, and `q` procedures use optional `log_p` and `lower_tail`
arguments where the upstream API does. Random procedures return allocatable
arrays. Multivariate random procedures return arrays whose rows are draws.

## Build

With FPM:

```text
fpm test
fpm run --example demo_extra_distr
```

With GNU Make:

```text
make check
make release
make MODE=release example
```

The checked build enables bounds, allocation, and other runtime checks. The
release build uses `-O3`. Both builds treat warnings as errors.

## Scope and differences from R

The numerical distributions are translated. R vector recycling, labels for
categorical values, data-frame/matrix coercion, warning collection, and Rcpp
registration are not reproduced. Fortran uses explicit scalar arguments and
type-safe arrays.

Three upstream edge behaviors are retained by default and corrected behavior
can be selected explicitly:

- `dslash(..., source_compatible=.false.)` applies the missing scale factor at
  the exact center of a non-unit-scale slash distribution.
- `pdlaplace(..., source_compatible=.false.)` branches relative to `location`
  rather than relative to zero.
- `ddirichlet(..., source_compatible=.false.)` enforces that the input lies on
  the probability simplex.

See `PORTING_NOTES.md`, `API_MAP.md`, and `VALIDATION.md` for details.

## License

GPL-2.0-only, matching the upstream `License: GPL-2` declaration. See
`LICENSE` and `NOTICE`.
