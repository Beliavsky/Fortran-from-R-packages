# optimflex-fortran

Modern Fortran/FPM translation of the computational code in the R package
`optimflex` 0.1.8.

## Included algorithms

- damped BFGS with strong-Wolfe line search
- limited-memory BFGS with box constraints (`l_bfgs_b`)
- pure Newton-Raphson
- modified Newton with dynamic diagonal ridge rescue
- Gauss-Newton
- Levenberg-Marquardt
- Powell dogleg
- double dogleg
- forward, central, and Richardson numerical gradient/Hessian/Jacobian routines
- Cholesky positive-definiteness checks
- the eight AND-combined convergence criteria exposed by the R package
- optional observed Hessian and Gauss-Newton-curvature callbacks

The package is standalone and has no R, BLAS, LAPACK, or `numDeriv` runtime
dependency.

## Build

```text
fpm build
fpm test
```

Examples:

```text
fpm run --example rosenbrock_suite
fpm run --example nonlinear_least_squares
```

## Callback style

```fortran
use optimflex

type(optim_result) :: res
real(dp) :: x0(2)

x0 = [-1.2_dp, 1.0_dp]
call bfgs(x0, objective, res)
```

Analytic gradient and Hessian callbacks are optional. When absent, the selected
finite-difference method is used.

## License

MIT, matching the upstream package. The supplied upstream tree is retained
under `original/optimflex-master/` for provenance.
