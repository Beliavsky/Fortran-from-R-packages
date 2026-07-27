# DEoptimR Modern Fortran

A modern Fortran 2018 translation of the computational algorithms in the R
package **DEoptimR 1.2-0**. The project provides callback-based global
optimization without R objects or runtime dependencies.

## Implemented numerical algorithms

- `jde_optimize`: sequential/asynchronous self-adaptive jDE
- `spjde_optimize`: synchronous population jDE corresponding to the numerical
  algorithm used by `SPJDEoptim`
- `ncde_optimize`: neighborhood-based niching differential evolution for
  multimodal optimization
- DE/rand/1/either-or/bin mutation and recombination
- self-adaptation of `F`, `CR`, mutation-choice probability, and NCDE
  neighborhood size
- dither and optional componentwise jitter
- deterministic midpoint/bounce-back bound handling
- inequality constraints and equality constraints converted with scalar or
  per-equality tolerances
- adaptive total-constraint-violation threshold
- injected initial candidate columns
- median- or maximum-based stopping criteria
- saved final populations, costs, constraints, and violations
- multimodal solution archives, fixed or automatically identified niche radius,
  archive replacement, and optional population reinitialization
- reproducible seeding through `seed_rng`

The callback interfaces are:

```fortran
function objective_function(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
end function

subroutine constraint_function(x, values)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: values(:)
end subroutine
```

See `example/constrained_example.f90` for a complete constrained call.

## Build and test

```sh
make debug
make release
make check
```

The debug configuration enables bounds and runtime checking. The release
configuration uses `-O2`. Both treat compiler warnings as errors.

An `fpm.toml` manifest is included. Validation for this release used the
provided Makefile/scripts because `fpm` was not installed in the validation
environment.

## Differences from the R package

- R functions and arbitrary `...` data are represented by typed Fortran
  procedure callbacks. Callback state can be held in a module.
- `spjde_optimize` preserves synchronous generation and selection semantics,
  but callback evaluation is serial. The R package's `mirai` process
  orchestration is not reproduced.
- Random streams use Fortran `random_number`; exact R random-stream and
  iteration-by-iteration equivalence is not claimed.
- NCDE archive-neighbor reinitialization is persisted into the next population,
  matching the documented intent. The R source updates the current population
  during construction of a separate next population, which can discard that
  reinitialization at the generation boundary.
- The R implementation's apparent use of the loop target index instead of the
  selected nearest-neighbor index in one constrained NCDE branch is corrected.
- Console formatting, R argument introspection, and R list/S3 infrastructure
  are excluded.

## License

The original package declares `GPL (>= 2)`. This translation is licensed under
**GPL-2.0-or-later**. Every Fortran source carries an SPDX identifier and GPL
notice. See `LICENSE` and `ORIGIN.md`.
