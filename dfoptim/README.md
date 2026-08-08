# dfoptim-fortran

A modern Fortran 2018 translation of the computational algorithms in the R
package `dfoptim` 2023.1.0.

## Implemented algorithms

- `hjk`: unbounded Hooke-Jeeves pattern search
- `hjkb`: box-constrained Hooke-Jeeves search
- `nmk`: modified Nelder-Mead with stagnation restarts
- `nmkb`: transformed-bound modified Nelder-Mead
- `mads`: randomized lower-triangular mesh-adaptive direct search

All optimizers are derivative-free. Objective functions are native Fortran
procedure callbacks and may receive optional polymorphic user data.

## Build with FPM

```text
fpm build
fpm test
fpm run --example dfoptim_example
```

## Minimal example

```fortran
program example
   use dfoptim, only : dp, nmk, nmk_control_t, dfoptim_result_t
   implicit none

   type(nmk_control_t) :: control
   type(dfoptim_result_t) :: result

   control%maxfeval = 10000
   result = nmk([-1.2_dp, 1.0_dp], rosenbrock, control)

   print '(a,2f14.8)', 'x = ', result%x
   print '(a,es14.6)', 'f = ', result%value

contains

   function rosenbrock(x, user_data) result(value)
      real(dp), intent(in) :: x(:)
      class(*), intent(inout), optional :: user_data
      real(dp) :: value

      value = 100.0_dp * (x(2) - x(1)**2)**2 + (1.0_dp - x(1))**2
   end function rosenbrock
end program example
```

## Public types

- `hj_control_t`
- `nmk_control_t`
- `mads_control_t`
- `dfoptim_result_t`
- `mads_log_t`

`dfoptim_result_t` contains the final parameter vector, objective value,
convergence code, evaluation and iteration counts, restart count, a message,
and the MADS iteration log where applicable.

## Bounds

`hjkb` and `nmkb` take elementwise lower and upper arrays. `nmkb` supports
finite, lower-only, upper-only, and unbounded components. A starting value must
be strictly inside every bound that is represented through a logarithmic or
hyperbolic transformation.

The current `mads` algorithm, like the supplied R implementation, accepts
either fully finite box bounds or a completely unbounded problem. Partially
finite MADS bounds are rejected.

## Reproducibility

Hooke-Jeeves direction ordering and MADS polling use a package-local random
number generator. Set `control%seed` to obtain reproducible runs without
modifying Fortran's global `random_number` state.

## Licensing

GPL-2.0-or-later. See `COPYING` and `NOTICE.md`.
