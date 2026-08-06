# lbfgs-fortran

Modern Fortran 2018 translation of the computational core of the R package
`lbfgs` 1.2.1.2.

The package provides native limited-memory BFGS optimization and orthant-wise
limited-memory quasi-Newton optimization. It does not require R, Rcpp, C, or
C++ at build or run time.

## Implemented algorithms

- L-BFGS two-loop inverse-Hessian recursion
- More-Thuente line search
- Backtracking Armijo line search
- Backtracking regular-Wolfe line search
- Backtracking strong-Wolfe line search
- OWL-QN for L1-regularized objectives
- Optional partial L1 regularization over a contiguous variable range
- Gradient-norm and past-objective stopping tests
- Optional per-iteration progress and cancellation callback
- Combined objective/gradient callbacks
- Separate objective and gradient callbacks
- Optional polymorphic user data
- Original libLBFGS status values and status messages

## Build with FPM

```text
fpm build
fpm test
fpm run --example rosenbrock
fpm run --example owlqn_soft_threshold
```

## Basic use

```fortran
program example
    use lbfgs
    implicit none

    real(dp) :: x(2)
    type(lbfgs_parameter_t) :: parameters
    type(lbfgs_result_t) :: result

    x = [-1.2_dp, 1.0_dp]
    parameters = lbfgs_parameter_t()
    parameters%epsilon = 1.0e-9_dp

    call lbfgs_minimize(evaluate, x, result, parameters)
    print *, x, result%value, result%message

contains

    subroutine evaluate(x, f, g, step, user_data)
        real(dp), intent(in) :: x(:), step
        real(dp), intent(out) :: f, g(:)
        class(*), intent(inout), optional :: user_data

        f = 100.0_dp * (x(2) - x(1)**2)**2 + (1.0_dp - x(1))**2
        g(1) = -400.0_dp * x(1) * (x(2) - x(1)**2) - &
            2.0_dp * (1.0_dp - x(1))
        g(2) = 200.0_dp * (x(2) - x(1)**2)
    end subroutine evaluate
end program example
```

The `step` callback argument is the current line-search step. A callback may
ignore it. The optional `user_data` argument can hold a derived-type problem
object, arrays, or another application-defined value.

## OWL-QN indexing

The R and C APIs use a zero-based half-open interval for
`orthantwise_start:orthantwise_end`. The Fortran API uses the normal Fortran
convention: a one-based inclusive interval.

For example, to regularize variables 2 through 10:

```fortran
parameters%orthantwise_start = 2
parameters%orthantwise_end = 10
```

An `orthantwise_end` value of zero means the final variable.

## License

The translated package is GPL-2.0-or-later, matching the R package. The
libLBFGS-derived numerical kernel retains the original MIT copyright and
permission notice in `licenses/libLBFGS-MIT.txt`.

Original source files are retained under `upstream/lbfgs-master` for license,
attribution, and validation traceability.
