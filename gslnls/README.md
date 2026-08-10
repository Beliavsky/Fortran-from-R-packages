# gslnls-fortran

Modern Fortran/FPM translation of the computational core of the R package
`gslnls` 1.4.2.

The original package provides an R interface to GNU Scientific Library (GSL)
nonlinear least-squares routines.  This translation is standalone: it does not
require R, GSL, Matrix, BLAS, or LAPACK.

## Implemented numerical functionality

- Nonlinear least squares with typed model callbacks.
- Analytic or finite-difference Jacobians.
- Levenberg-Marquardt (`NLS_LM`).
- Levenberg-Marquardt with geodesic acceleration (`NLS_LMACCEL`).
- Powell dogleg (`NLS_DOGLEG`).
- Double dogleg (`NLS_DDOGLEG`).
- Two-dimensional subspace trust-region step (`NLS_SUBSPACE2D`).
- Steihaug-Toint truncated conjugate gradient (`NLS_CGST`).
- Matrix-free large-system Steihaug-CG through a Jacobian/Jacobian-transpose
  operator callback.
- Parameter lower and upper bounds.
- Vector weights and general positive-definite weight matrices.
- Robust IRLS regression with the losses exposed by `gsl_nls_loss`:
  Huber, Barron, bisquare, Welsh, optimal, Hampel, GGW, and LQQ.
- Multi-start optimization with quasi-random sampling and local refinement.
- Forward and central finite-difference Jacobians.
- Directional second derivatives for geodesic acceleration.
- Covariance, hat values, Cook distances, Gaussian log likelihood, and normal
  confidence intervals.
- Iteration traces and evaluation counters.

## Quick example

```fortran
program demo
  use gslnls, only : dp, nls_result, fit_nls
  implicit none
  real(dp), parameter :: x(5) = [0._dp,1._dp,2._dp,3._dp,4._dp]
  real(dp), parameter :: y(5) = [1._dp,3._dp,5._dp,7._dp,9._dp]
  type(nls_result) :: fit

  call fit_nls(model, y, [0._dp,0._dp], fit, jac=jacobian)
  print *, fit%par
contains
  subroutine model(par,yhat,ierr)
    real(dp), intent(in) :: par(:)
    real(dp), intent(out) :: yhat(:)
    integer, intent(out) :: ierr
    yhat = par(1) + par(2)*x
    ierr = 0
  end subroutine model

  subroutine jacobian(par,j,ierr)
    real(dp), intent(in) :: par(:)
    real(dp), intent(out) :: j(:,:)
    integer, intent(out) :: ierr
    j(:,1) = 1._dp
    j(:,2) = x
    ierr = 0
    if (size(par) /= 2) ierr = 1
  end subroutine jacobian
end program demo
```

## Build

```text
fpm build
fpm test
```

Strict GNU Fortran checks can be run with:

```text
scripts/test_gfortran.sh
```

or on Windows:

```text
scripts\test_gfortran.bat
```

The strict scripts use:

```text
-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

## Design notes

The public API works directly with arrays and procedure callbacks instead of R
formulas, environments, sparse Matrix objects, or S3 classes. User callback
invocations in the library are made from module procedures with explicit
procedure interfaces. This is intentional for portability with gfortran builds
that promote implicit-interface warnings to errors.

The package-owned multi-start logic is retained, but this first Fortran release
uses a Halton sequence for every dimension. Upstream uses GSL Sobol points below
41 parameters and Halton points above that threshold. See
`TRANSLATION_COVERAGE.md` for all fidelity notes.

## Licensing

The upstream package declares `LGPL-3`. Translated source files are marked
`LGPL-3.0-only`, and the GNU LGPL v3 text is included in `LICENSE`. The complete
supplied upstream package is retained under `original/gslnls-master/`.
