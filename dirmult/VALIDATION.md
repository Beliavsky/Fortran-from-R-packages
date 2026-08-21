# Validation

Validation was performed with GNU Fortran 14.2.0 using:

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -fbacktrace
```

FPM was not installed in the execution environment, so the exact source and
test units described by `fpm.toml` were compiled and linked manually with
`gfortran`.

## Tests

- `test_core`: Weir-Hill MoM, observed/expected Fisher information, observed-
  and expected-information MLE paths, zero-row/zero-column filtering,
  transformed standard errors, and multinomial log-likelihood.
- `test_profiles`: fixed-theta profile likelihood, profile grid, and a two-table
  common-theta constrained fit.
- `test_simulation`: Dirichlet simulation moments, Dirichlet-multinomial row
  totals, and the simulation null-test result path.

The core MLE reference was independently checked with SciPy optimization of the
same Dirichlet-multinomial log-likelihood. For the test table:

- gamma = approximately `(1.03582212, 1.15481351, 0.98728973)`
- pi = approximately `(0.32594287, 0.36338598, 0.31067115)`
- theta = `0.2393532469974`
- log-likelihood kernel = `-227.0699809395`

The common-theta test solves two independent tables jointly and converges to
approximately `theta = 0.1942962982040` with combined log-likelihood kernel
`-464.5865698622`.

All tests pass with array bounds, allocation, and other runtime checks enabled.
All Fortran source/test/example lines are at most 132 characters, so no
compiler-specific free-line-length option is required.
