# roptim-fortran

Version 0.1.1.

A modern Fortran 2018 computational port of the R package **roptim 0.1.7**.
The package provides one callback-based interface to five optimization methods:

- Nelder-Mead
- BFGS
- nonlinear conjugate gradients
- L-BFGS-B
- simulated annealing (SANN)

The original R package is a C++ wrapper around the optimization routines used
by `stats::optim()`. This port removes R, Rcpp, and Armadillo and supplies
native Fortran implementations and typed result/control objects.

## Build with FPM

```console
fpm build
fpm test
fpm run --example rosenbrock_methods
fpm run --example wild_sann
```

## Minimal example

```fortran
program demo
  use roptim_mod, only : dp, roptim_control_t, roptim_result_t, &
       roptim_minimize, method_bfgs
  implicit none

  real(dp) :: x(2)
  type(roptim_control_t) :: control
  type(roptim_result_t) :: result

  x = [-1.2_dp, 1.0_dp]
  control%compute_hessian = .true.
  call roptim_minimize(x, objective, result, method_bfgs, &
       gradient=objective_gradient, control=control)

  print *, x, result%value, result%success

contains

  function objective(x, user_data) result(f)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: f
    f = 100.0_dp*(x(2)-x(1)**2)**2 + (1.0_dp-x(1))**2
  end function objective

  subroutine objective_gradient(x, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    class(*), intent(inout), optional :: user_data
    g(1) = -400.0_dp*x(1)*(x(2)-x(1)**2)-2.0_dp*(1.0_dp-x(1))
    g(2) = 200.0_dp*(x(2)-x(1)**2)
  end subroutine objective_gradient

end program demo
```

If no gradient procedure is supplied, central finite differences are used.
L-BFGS-B finite differences are bound-aware.

## Main API

- `roptim_minimize`
- `roptim_approximate_gradient`
- `roptim_approximate_hessian`
- `roptim_control_t`
- `roptim_result_t`

The control type follows the names used by R `optim()` where practical:
`fnscale`, `parscale`, `ndeps`, `reltol`, `abstol`, Nelder-Mead coefficients,
CG type, L-BFGS-B memory/`factr`/`pgtol`, and SANN temperature/`tmax`.

## Licensing

The port retains the original package's GPL-2-or-later license. The bundled
L-BFGS-B 3.0 numerical kernel retains its separately distributed BSD-style
license and notice. See `LICENSE`, `LICENSES/`, and `NOTICE.md`.

## Compatibility notes

The numerical API is native Fortran, not a source-compatible translation of
C++ templates. Floating-point paths and evaluation counts are not guaranteed
to be bit-identical to R's `optim()`, particularly for BFGS, CG, and
Nelder-Mead. See `docs/PORTING_NOTES.md` and `docs/API_MAP.md`.
