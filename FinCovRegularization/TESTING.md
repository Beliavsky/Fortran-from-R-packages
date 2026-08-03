# Testing

## FPM

```text
fpm test
```

## GNU Fortran scripts

On Unix-like systems:

```text
./run_gfortran_tests.sh
```

On Windows with GNU Fortran available on `PATH`:

```text
run_gfortran_tests.bat
```

The strict debug build uses Fortran 2018, warnings as errors, run-time bounds
and argument checking, and floating-point traps for invalid operations,
division by zero, and overflow.

The tests cover:

- Frobenius and operator norms
- All four regularization operators and `threshold_min`
- Independence covariance
- One- and multi-factor macro models
- Fundamental OLS and WLS models
- Statistical factor reconstruction and automatic factor selection
- Short-enabled, long-only, and boundary GMVP solutions
- Risk parity
- Banding, tapering, hard-threshold, and soft-threshold CV
- Deterministic repeated splits and operator-norm CV loss
