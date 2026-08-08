# mcga-fortran

Modern Fortran/FPM translation of the computational code in the R package
`mcga` 3.0.9 (Machine Coded Genetic Algorithms for Real-Valued Optimization
Problems).

The upstream package is licensed GPL-2.0-or-later. This translation is
released under the same terms. The complete input package is retained under
`original/mcga-master/` for provenance.

## Included

- Single-objective machine-coded genetic algorithm (`mcga_optimize`), matching
  the exported R `mcga()` engine.
- Multi-objective machine-coded genetic algorithm (`multi_mcga_optimize`),
  including the upstream rank-score calculation and tournament behavior.
- Native byte representation conversion for `real(real64)` values.
- One-point, two-point, and uniform byte crossover.
- Byte mutation by +/-1 with unsigned-byte wraparound.
- Random-byte mutation.
- Bound repair by uniform resampling.
- `byte_mutation`, dynamic byte mutation, random byte mutation, and their
  dynamic variants.
- Byte, one-point-byte, two-point-byte, SBX, flat, arithmetic, BLX, linear,
  and unfair-average crossover routines.
- Deterministic random seeding helper for reproducible Fortran tests/examples.

## Not included

`mcga2()` is a wrapper around the external R package `GA::ga`; the GA package's
selection/population/termination engine is not source code belonging to mcga
and is therefore not reimplemented here. All mcga-specific operators used by
`mcga2()` are translated and can be plugged into a Fortran GA driver.

R `.Call`/Rcpp registration, S4/GA objects, R environments, parallel R
execution, monitors, and documentation-only plotting are omitted.

See `TRANSLATION_COVERAGE.md` for exact source-fidelity notes.

## Build

```sh
fpm build
fpm test
```

Examples:

```sh
fpm run --example mcga_example
fpm run --example multi_mcga_example
```

The library has no external numerical dependencies.

## Basic use

```fortran
program demo
  use mcga, only : dp, mcga_result, mcga_optimize
  implicit none
  type(mcga_result) :: res

  call mcga_optimize(100, 2, -10.0_dp, 10.0_dp, objective, res, &
                     maxiter=200, seed=123)
  print *, res%population(1,:), res%costs(1)
contains
  function objective(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = sum((x - [2.0_dp, -1.0_dp])**2)
  end function objective
end program demo
```

## Byte representation

MCGA works on the in-memory bytes of IEEE-style floating-point values. The
Fortran implementation uses `TRANSFER` and 0--255 integer byte values, so the
algorithm follows the host machine's native `real(real64)` representation and
byte order, just as the C/C++ original follows the host C `double`
representation and byte order.
