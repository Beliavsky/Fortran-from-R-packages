# rangen-fortran

Modern Fortran/FPM translation of the computational code in the R package `rangen` 0.0.1.

The upstream package provides fast random-number generators, matrix/column generators, sampling utilities, seeding, and a high-resolution timer. This port removes the Rcpp/RcppArmadillo interface layer and exposes the same numerical ideas directly from Fortran.

## Implemented functionality

### Random generators

- `runif`
- `rbeta`
- `rexp`
- `rchisq`
- `rgamma`
- `rgeom`
- `rcauchy`
- `rt`
- `rpareto`
- `rfrechet`
- `rlaplace`
- `rgumbel` and the upstream-compatible misspelling `rgumble`
- `rarcsine`
- `rnorm`

Each generator has a scalar-parameter matrix form such as `rnorm_mat`, and the R `colR*` interfaces are represented by `col_r*` functions accepting parameter vectors with R-style recycling.

### Sampling and utilities

- `sample_int`
- `sample_real`
- `col_sample`
- `row_sample`
- `set_seed`
- `nano_time`

`set_seed` seeds the distribution, uniform, normal, and sampling streams together. The PCG32 recurrence from the upstream package is implemented with 16-bit limb arithmetic, avoiding reliance on signed 64-bit integer overflow.

## Example

```fortran
program demo
    use rangen
    implicit none
    real(dp), allocatable :: x(:)
    integer, allocatable :: ix(:)

    call set_seed(42_i8)
    x = rnorm(10000, 2.0_dp, 3.0_dp)
    print *, sum(x) / real(size(x), dp)

    ix = sample_int(10, 5, .false.)
    print *, ix
end program demo
```

## Build

With FPM:

```text
fpm test
fpm run --example basic
```

The project has no external runtime dependencies.

## Validation

The source tree is tested with GNU Fortran using:

```text
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all
```

The tests cover the PCG32 recurrence, seed reproducibility, statistical moments for all distribution families, matrix/column generation, and sampling with and without replacement.

## Important compatibility notes

The target distributions and public computational behavior are preserved, but the stream is not intended to be bit-for-bit identical to the R package for every generator. The upstream normal generator comes from the external `zigg` package, whose implementation is not bundled in the supplied source. This port uses a Box-Muller normal generator driven by an independent PCG32 stream.

Several clear upstream defects are corrected; see `PORTING_NOTES.md`.

## License

The upstream package declares `GPL-3`. This translation is distributed under GPL-3.0-only. See `LICENSE`, `LICENSES.md`, and the preserved upstream snapshot under `upstream/`.
