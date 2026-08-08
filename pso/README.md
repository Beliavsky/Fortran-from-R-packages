# pso-fortran

Modern Fortran/FPM translation of the computational portions of the R package
`pso` 1.0.4 by Claus Bendtsen.

The upstream package implements particle swarm optimization consistent with
Standard PSO 2007 and Standard PSO 2011.  This translation keeps the numerical
optimizer, benchmark functions, repeated benchmark runner, success-rate curve,
and efficiency/statistics calculations.  R S4 dispatch and plotting methods are
not translated.

## Features

- SPSO2007 and SPSO2011
- scalar or vector lower/upper bounds
- optional supplied initial point (use IEEE NaNs to request random initialization)
- random or fixed particle processing order
- informant topology and topology regeneration
- synchronous (`vectorize`) and asynchronous particle updates
- velocity clamping
- restart-on-small-swarm-diameter logic
- stagnation and function-evaluation stopping criteria
- minimization or maximization through `fnscale`
- optional hybrid bounded limited-memory BFGS refinement
- analytical or numerical local gradients
- trace-history storage
- upstream Parabola, Griewank, shifted Rosenbrock, Rastrigin, and Ackley tests
- repeated benchmark summaries and success/efficiency calculations

## Build

```text
fpm build
fpm test
```

Examples:

```text
fpm run --example rastrigin_example
fpm run --example hybrid_example
```

A direct GNU Fortran validation script is also supplied:

```text
scripts/test_gfortran.sh
```

On Windows:

```text
scripts\test_gfortran.bat
```

## Minimal example

```fortran
program demo
   use pso, only : dp, pso_control, pso_result, psoptim, seed_random
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   implicit none

   real(dp) :: par(2)
   type(pso_control) :: control
   type(pso_result) :: result

   call seed_random(1)
   par = ieee_value(0.0_dp, ieee_quiet_nan)
   control%abstol = 1.0e-8_dp

   call psoptim(par, objective, -5.0_dp, 5.0_dp, result, control)
   print *, result%par, result%value

contains

   function objective(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = sum(x*x)
   end function objective
end program demo
```

See `API.md` and `TRANSLATION_COVERAGE.md` for details.

## License

The upstream package declares `License: LGPL-3`.  The translated code is
provided under LGPL-3.0, with SPDX identifiers in translated source files.
`COPYING.LESSER` contains the GNU Lesser General Public License version 3.
The original package source is retained under `original/pso-master/` for
provenance.
