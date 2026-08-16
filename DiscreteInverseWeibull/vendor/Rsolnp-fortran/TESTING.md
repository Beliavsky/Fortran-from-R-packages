# Testing

## FPM

```text
fpm test
fpm run
fpm run --example basic_constrained
fpm run --example multistart_example
fpm run --example standard_form_example
```

## Direct GNU Fortran validation

On Unix-like systems:

```text
./scripts/test_gfortran.sh
```

On Windows with GNU Fortran available on `PATH`:

```text
scripts\test_gfortran.bat
```

The scripts compile modules in dependency order and run all tests, examples,
and the demo. Set `FC` to select a compiler executable.

## Test programs

- `test_solver`: analytic and finite-difference Rosenbrock optimization.
- `test_constraints`: equality, inequality, and box-constrained benchmarks.
- `test_multistart`: deterministic starts, bound compliance, and HS05 multistart.
- `test_standardize_benchmarks`: registry coverage, standard-form expansion,
  Powell solution, and KKT feasibility.

Validation is intended for both a runtime-checked configuration and an
optimized warning-as-error configuration. Floating-point traps are enabled in
the strict GNU build.
