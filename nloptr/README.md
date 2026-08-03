# nloptr-fortran

A self-contained modern Fortran/FPM translation of the computational interface
of the R package `nloptr` 2.2.1.9000.

The library solves problems of the form

```text
minimize f(x)
subject to g(x) <= 0
           h(x)  = 0
           lower <= x <= upper
```

Objective and constraint functions are supplied as typed Fortran procedure
callbacks. Results are returned in a `nloptr_result` derived type.

## Implemented public API

All 26 exported R entry points have Fortran counterparts:

- Main interface and options: `nloptr`, `nl_opts`,
  `nloptr_get_default_options`, `nloptr_print_options`, `is_nloptr`
- Derivatives: `nl_grad`, `nl_jacobian`, `check_derivatives`
- Local gradient methods: `lbfgs`, `varmetric`, `tnewton`, `slsqp`, `mma`,
  `ccsaq`, `auglag`
- Local derivative-free methods: `neldermead`, `sbplx`, `cobyla`, `bobyqa`,
  `newuoa`
- Global/multistart methods: `direct`, `direct_l`, `crs2lm`, `isres`, `stogo`,
  `mlsl`

Dots in R names are replaced by underscores where required by Fortran.

## Numerical implementation

This is a translation of the high-level package interface, not a line-by-line
rewrite of the bundled NLopt 2.10.0 C/C++ library. The self-contained backend
uses:

- projected full-memory BFGS with safeguarded Armijo line search;
- Nelder-Mead simplex search;
- coordinate trust-region/pattern search;
- quadratic-penalty handling of equality and inequality constraints;
- deterministic Halton multistart for global wrappers;
- central finite-difference gradients and Jacobians.

The algorithm-name mapping is documented in `PORTING.md`. Applications needing
bit-for-bit NLopt behavior or its specialized asymptotic performance should link
to the native NLopt library instead.

## Minimal example

```fortran
program example
  use nloptr_mod
  use nloptr_example_functions, only: rosenbrock_objective
  implicit none

  type(nloptr_problem) :: problem
  type(nloptr_options) :: options
  type(nloptr_result) :: result

  problem%n = 2
  problem%objective => rosenbrock_objective
  allocate(problem%lower(2), problem%upper(2))
  problem%lower = -huge(1.0_dp)
  problem%upper =  huge(1.0_dp)

  options = nl_opts('NLOPT_LD_LBFGS')
  options%maxeval = 5000
  call nloptr(problem, [-1.2_dp, 1.0_dp], options, result)

  print *, result%solution
  print *, result%objective
end program example
```

## Build

```text
fpm build
fpm test
fpm run --example local_optimization
fpm run demo_nloptr
```

Direct GNU Fortran validation scripts are also provided in `scripts/`.
