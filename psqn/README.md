# psqn-fortran

Modern Fortran/FPM translation of the computational code in the R package
**psqn 0.3.2** (Partially Separable Quasi-Newton), originally by Benjamin
Christoffersen.

The port preserves the upstream Apache license, removes R/Rcpp/Eigen plumbing,
and skips non-computational presentation/plotting material. It is intended as a
standalone numerical library for Fortran programs.

## Implemented algorithms

* Partially separable quasi-Newton optimization with arbitrary element maps.
* Structured global/private parameter layout corresponding to the main psqn API.
* Per-element BFGS or SR1 Hessian updates.
* Newton directions by preconditioned conjugate gradients.
* No, diagonal, regularized-Cholesky, and structured block preconditioning.
* Weak or strong Wolfe line search with interpolation/zoom logic.
* Parameter masking.
* Full inverse-BFGS optimizer.
* Equality constraints with an augmented-Lagrangian outer iteration.
* Private-block optimization for structured problems.
* Numerical Hessian assembly using Richardson extrapolation.
* Packed symmetric linear-algebra helpers and compensated summation.

See `TRANSLATION_NOTES.md` for exact upstream-to-Fortran source mapping and
intentional implementation differences.

## Build and test with FPM

```text
fpm build
fpm test
```

Examples can be run with

```text
fpm run --example generic_quadratic
fpm run --example structured_quadratic
```

The source also builds without external libraries using a Fortran 2018 compiler.
During translation it was checked with GNU Fortran using strict interface and
runtime checking, and all included tests passed.

## Generic partially separable example

```fortran
program generic_quadratic
  use psqn
  implicit none

  type(psqn_element_spec), allocatable :: specs(:)
  type(psqn_info) :: info
  real(dp) :: x(3)

  allocate(specs(2))
  specs(1)%idx = [1, 2]
  specs(2)%idx = [2, 3]
  x = 0.0_dp

  call psqn_optimize_generic(x, specs, element, info)
  print '(a,3f12.6)', 'x = ', x
  print '(a,es14.6)', 'f = ', info%value

contains

  subroutine element(i, z, f, g, comp_grad)
    integer, intent(in) :: i
    real(dp), intent(in) :: z(:)
    real(dp), intent(out) :: f
    real(dp), intent(out) :: g(:)
    logical, intent(in) :: comp_grad

    select case (i)
    case (1)
      f = 0.5_dp * ((z(1) - 1.0_dp)**2 + (z(2) + 2.0_dp)**2)
      if (comp_grad) g = [z(1) - 1.0_dp, z(2) + 2.0_dp]
    case (2)
      f = 0.5_dp * ((z(1) + 2.0_dp)**2 + (z(2) - 3.0_dp)**2)
      if (comp_grad) g = [z(1) + 2.0_dp, z(2) - 3.0_dp]
    end select
  end subroutine element

end program generic_quadratic
```

The minimizer is `(1, -2, 3)`.

## Structured layout

For `global_dim = G` and private dimensions `p(1:m)`, element `i` receives

```text
[ global parameters 1:G, private parameters belonging to element i ]
```

Use `psqn_optimize_structured` when this is the natural structure. Use
`psqn_make_structured_specs` if explicit generic element maps are useful.

## Licensing

See `LICENSE`, `NOTICE`, and `original/DESCRIPTION`. The translated source files
carry `SPDX-License-Identifier: Apache-2.0`.
