# onls-fortran

Modern Fortran/FPM translation of the computational core of R package
`onls` 0.1-4 (Orthogonal Nonlinear Least-Squares Regression).

## What is translated

The package implements the two-stage numerical algorithm used by upstream
`onls()`:

1. ordinary vertical nonlinear least squares to obtain starting parameters;
2. orthogonal nonlinear least squares, where each residual is the Euclidean
   distance from an observed `(x,y)` point to its nearest point on the model
   curve.

The outer nonlinear least-squares problems use a standalone bounded
Levenberg-Marquardt/Gauss-Newton solver.  Each nearest-point problem uses Brent
one-dimensional minimization over the same global/local intervals used by the
R source.

Implemented features include:

- scalar predictor / scalar response nonlinear models via an explicit callback;
- ordinary-NLS initialization followed by orthogonal-NLS fitting;
- lower/upper parameter bounds;
- observation weights;
- fixed parameters in the orthogonal stage;
- `window` and `extend` projection intervals;
- `x0`/`y0` nearest points;
- vertical and orthogonal residuals/deviances;
- numerical parameter-gradient matrix, covariance, and standard errors;
- vertical and orthogonal Gaussian log-likelihood helpers;
- the upstream orthogonality-angle diagnostic;
- NIST Chwirut2 example data.

The library is standalone and requires no R, `minpack.lm`, BLAS, or LAPACK.

## Basic use

```fortran
program demo
    use onls
    implicit none
    real(dp), parameter :: x(6) = [-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp,3.0_dp]
    real(dp), parameter :: y(6) = [-4.2_dp,-0.7_dp,2.3_dp,5.4_dp,7.7_dp,11.2_dp]
    type(onls_result) :: fit

    call fit_onls(line_model, x, y, [0.0_dp,1.0_dp], fit)
    print *, fit%par_onls
contains
    subroutine line_model(xx, par, yy, ierr)
        real(dp), intent(in) :: xx(:), par(:)
        real(dp), intent(out) :: yy(:)
        integer, intent(out) :: ierr
        yy = par(1) + par(2)*xx
        ierr = 0
    end subroutine line_model
end program demo
```

## Build

```text
fpm build
fpm test
fpm run --example line_example
fpm run --example chwirut2_example
```

For strict GNU Fortran validation without FPM, use
`scripts/test_gfortran.sh` or `scripts/test_gfortran.bat`.

## Important fidelity notes

The R package delegates both nonlinear fits to `minpack.lm::nls.lm`; this
translation supplies its own bounded LM implementation, so iteration counts and
last-bit parameter values are not expected to be identical to MINPACK.
The orthogonal objective and projection intervals follow the R source.

The upstream source sorts predictor/response values before fitting but does not
reorder the already-created weight vector.  `onls_control%mimic_r_unsorted_weights`
is `.true.` by default to reproduce that behavior.  Set it to `.false.` for the
more conventional interpretation in which weights follow their observations
through sorting.

See `TRANSLATION_COVERAGE.md` for the full mapping and omissions.
