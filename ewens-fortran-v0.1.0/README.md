# ewens-fortran

Modern free-format Fortran translation of the computational code in the R
package **ewens 0.1.0**.

The package implements the Ewens sampling distribution, Chinese restaurant
process sampling, the generalized Pitman-Yor/Chinese restaurant process,
GEM stick breaking, the distribution of the number of classes, and maximum
likelihood estimation of the diversity parameter.

## Features

- `dewens()` and `dewens_log()` for class-membership vectors.
- `dewens_counts()` and `dewens_counts_log()` for frequency spectra
  `m(j) = number of classes observed j times`.
- `dewens_k()` / `dewens_k_log()` for the number-of-classes distribution.
- `ewens_k_exact()` for the exact expected number of classes.
- `rewens()` for Ewens/Chinese-restaurant sampling.
- `gcrp()` for the generalized Chinese restaurant (Pitman-Yor) process.
- `rgem()` for GEM stick-breaking weights.
- `ewens_mle()` and `ewens_mle_nk()` for maximum-likelihood estimation.
- `ewens_seed()` for reproducible simulations.

The attached `copula-fortran` translation is used as a local FPM dependency.
Its exact `stirling_first()` routine is used for `dewens_k()` through `n=20`.
For larger `n`, `ewens-fortran` switches to an overflow-safe log-space
Stirling recurrence.

## Build

```text
fpm build
fpm test
fpm run --example ewens_demo
```

All compiled Fortran source is free-format `.f90`. No C code is compiled or
called. The original C/R sources under `upstream/` and the dependency's own
provenance directories are retained only for attribution and auditability.

## Example

```fortran
program example
  use ewens, only : dp, i8, ewens_seed, rewens, dewens, ewens_mle
  implicit none
  integer, allocatable :: x(:)

  call ewens_seed(2923_i8)
  x = rewens(24, 1.0_dp)
  print *, dewens(x, 1.0_dp)
  print *, ewens_mle(x)
end program example
```

## License

The `ewens`-derived translation source is MIT-licensed, preserving upstream
terms. The supplied `copula-fortran` dependency is GPL-3.0-or-later; therefore
linked combined binaries are subject to GPL-3.0-or-later-compatible terms.
See `LICENSE.md` and the dependency's own license files.
