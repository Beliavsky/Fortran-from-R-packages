# Testing

## FPM

```text
fpm test
fpm run
```

## Direct GNU Fortran validation

On Unix-like systems:

```text
./run_gfortran_tests.sh
```

On Windows:

```text
run_gfortran_tests.bat
```

The strict configuration uses Fortran 2018 conformance, warnings as errors,
runtime bounds checking, floating-point traps, and backtraces. Release
validation uses `-O3 -Werror`.

## Test programs

- `test_distributions`: densities, weighted statistics, matrix square roots,
  ellipses, and permutations
- `test_parametric`: normal, gamma, multivariate-normal, multinomial, and
  repeated-measures mixtures
- `test_regression`: linear, logistic, Poisson, expert-gated, and grouped
  regression mixtures
- `test_semiparametric_reliability`: product-kernel EM, symmetric mixtures,
  censored reliability models, model selection, bootstrap, and MCMC

Tests use deterministic synthetic data and the package's portable RNG.
