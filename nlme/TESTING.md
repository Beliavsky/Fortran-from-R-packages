# Testing

The test suite contains four programs:

- `test_covariance`: correlation, variance, ARMA-ACF, and positive-definite
  parameterizations.
- `test_gls`: REML GLS with estimated AR(1) correlation.
- `test_lme`: random-intercept LME, BLUPs, and conditional residuals.
- `test_nonlinear_diagnostics`: GNLS, nonlinear mixed-effects linearization,
  residual ACF, variogram, pooled SD, and simulation.

Run with FPM:

```text
fpm test
```

Direct GNU Fortran scripts are also included:

```text
./run_gfortran_tests.sh debug
./run_gfortran_tests.sh release
```

The debug configuration enables bounds checking, backtraces, and floating-point
traps. The release configuration uses `-O3`. Both treat warnings as errors.
