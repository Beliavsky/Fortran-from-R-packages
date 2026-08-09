# Testing

Run with FPM:

```sh
fpm test
```

Or with GNU Fortran:

```sh
scripts/test_gfortran.sh
```

The four tests cover:

1. Exact regression values from upstream `test1.R` and `test2.R`.
2. Native compact constraints and pre-factorized Hessian input.
3. Talbot-Katz equality/inequality cases for sizes 1 through 10.
4. Inconsistent constraints, non-positive-definite Hessians, invalid inputs,
   unconstrained problems, primal feasibility, and KKT stationarity.

The validation script compiles with Fortran 2018 conformance, warnings as
errors, bounds checking, uninitialized-use checks, and floating-point traps.
An optimized `-O3` configuration is also provided.
