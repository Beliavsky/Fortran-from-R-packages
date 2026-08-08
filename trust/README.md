# trust-fortran

Modern Fortran translation of the computational core of Charles J. Geyer's R package `trust` 0.1-9.

The package implements second-derivative trust-region local optimization using the same eigendecomposition-based trust-region subproblem described in the upstream R source.  It supports minimization, maximization, parameter rescaling, restricted objective domains represented by signed infinities, and optional iteration history.

## Build

```text
fpm build
fpm test
```

The library is standalone and does not require BLAS or LAPACK.  A self-contained cyclic Jacobi eigensolver is used for real symmetric Hessians, and a safeguarded bisection solve is used for the one-dimensional secular equation.

## Minimal use

```fortran
use trust

type(trust_options) :: options
type(trust_result) :: result

options%rinit = 1.0_dp
options%rmax = 5.0_dp
call trust_optimize(objective, x0, options, result)
```

The objective callback has the interface

```fortran
subroutine objective(x, value, gradient, hessian, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    real(dp), intent(out) :: gradient(:)
    real(dp), intent(out) :: hessian(:, :)
    integer, intent(out) :: status
end subroutine
```

Set `status=0` on a successful evaluation.  For a point outside a restricted domain, return positive infinity for minimization or negative infinity for maximization.  Gradient and Hessian values are ignored for such infeasible points.

See `API.md`, `TRANSLATION_COVERAGE.md`, and the programs under `example/`.

## License

The upstream package is MIT licensed (`MIT + file LICENSE` in DESCRIPTION).  Its full permission notice is retained in `LICENSE`, and the complete supplied package is retained under `original/trust-master/`.
