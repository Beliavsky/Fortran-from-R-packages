# contfrac-fortran

Modern Fortran/FPM translation of the computational code in the R package
`contfrac` 1.1-12 by Robin K. S. Hankin.

The upstream package evaluates ordinary and generalized continued fractions,
constructs their partial convergents, and expands a real number as a simple
continued fraction.  The original package implements the main recurrences in C
and exposes them through R wrappers.  This translation implements the same
algorithms directly in modern Fortran and has no external numerical dependency.

## Implemented computational API

- `cf`: ordinary continued fractions, real and complex.
- `gcf`: generalized continued fractions, real and complex, using the modified
  Lentz algorithm used by the upstream C implementation.
- `convergents`: numerator/denominator sequences for ordinary continued
  fractions, real and complex.
- `gconvergents`: numerator/denominator sequences for generalized continued
  fractions, real and complex.
- `as_cf`: simple continued-fraction expansion of a real scalar.
- `contfrac_info`: optional convergence metadata for `cf`/`gcf`.

Fortran is case-insensitive, so these routines correspond directly to the R
names `CF`, `GCF`, `convergents`, `gconvergents`, and `as_cf`.

## Build with FPM

```text
fpm build
fpm test
fpm run --example continued_fraction_demo
```

The package uses only standard Fortran 2018 facilities.

## Example

```fortran
program demo
    use contfrac, only : dp, cf, as_cf
    implicit none

    real(dp), allocatable :: a(:), terms(:)
    integer :: i

    a = [(1.0_dp, i = 1, 100)]
    print *, cf(a)                 ! golden ratio

    terms = as_cf(acos(-1.0_dp), 10)
    print *, terms                ! 3, 7, 15, 1, 292, ...
end program demo
```

For a finite generalized continued fraction, pass `finite=.true.`.  This mirrors
upstream behavior: the modified Lentz iteration may still terminate early when
machine-precision convergence is reached, while `finite=.true.` means that a
result is accepted even if the supplied finite sequence ends before the
convergence criterion is reached.

## Validation

The tests include direct translations of the upstream `tests/aaa.R` examples:
continued fractions for `sqrt(11)`, `sqrt(71)`, `exp(1)`, Euler's constant, and
the complex continued fraction for `tan(1+i)`.  Additional tests cover real and
complex partial convergents, infinite-coefficient truncation, nonconvergence
status, and rational `as_cf` termination.

## License and provenance

The upstream `DESCRIPTION` declares `License: GPL-2`.  This translation retains
that license.  See `COPYING` for the GNU GPL version 2 text and `UPSTREAM.md` for
provenance details.  Relevant original R/C source files are retained under
`upstream/` for auditability and attribution.
