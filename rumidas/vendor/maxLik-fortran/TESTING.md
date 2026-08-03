# Testing

## FPM

```text
fpm build
fpm test
fpm run
fpm run --example normal_mle
```

## GNU Fortran scripts

On Unix-like systems:

```text
./scripts/test_gfortran.sh
./scripts/test_gfortran_optimized.sh
```

On Windows with GNU Fortran available in `PATH`:

```text
scripts\test_gfortran.bat
```

The strict script enables warnings as errors, runtime array/bounds checks,
floating-point traps, and backtraces. The optimized script recompiles at `-O3`.

## Regression programs

- `test_derivatives`: analytic versus numerical derivatives.
- `test_deterministic_optimizers`: NR, BFGS, BFGSR, CG, and Nelder-Mead.
- `test_bhhh_and_inference`: normal MLE, observation scores, and robust covariance.
- `test_constraints_and_fixed`: equality constraints and fixed parameters.
- `test_stochastic`: deterministic seeded Adam and SGA convergence.
- `test_numeric_and_utilities`: objective-only numerical Newton, vector
  Jacobians, active-parameter utilities, and condition numbers.

The expected optima are closed-form or independently calculated. Tests do not
merely assert that routines return without error.
