# Validation

The translation was validated from a clean build directory with GNU Fortran
14.2.0 using:

```text
-std=f2018 -fcheck=all -Wall -Wextra -Werror=implicit-interface
```

FPM was not installed in the translation environment, so the exact source,
test, and example units described by `fpm.toml` were compiled and linked
manually.

Clean-build results:

```text
test_distributions: PASS
test_mcnode: PASS
test_sampling_dependence: PASS, rho=  0.800479
test_tornado_multivariate: PASS
test_workflow: PASS
```

Example output:

```text
total dimensions: 100 50 1
total type: VU
median of conditional means:   0.607000
```

Additional checks:

- no implicit-interface compilation errors;
- no source, test, or example lines longer than 132 characters;
- the supplied `mvtnorm-fortran` subtree is unchanged from the attached
  dependency archive;
- runtime array/bounds checking was enabled for all tests.

`-Wall -Wextra` reports expected warnings for exact floating-point comparisons
used at distribution support boundaries, discrete/tie comparisons, and similar
test assertions. These do not prevent compilation or test completion.
