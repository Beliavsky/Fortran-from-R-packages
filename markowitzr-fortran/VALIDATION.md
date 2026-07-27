# Validation

## Environment

- GNU Fortran 14.2.0
- Fortran 2018 mode
- Linux x86-64

Strict debug flags:

```text
-std=f2018
-Wall
-Wextra
-Wpedantic
-Wconversion-extra
-Wimplicit-interface
-Werror
-fcheck=all
-fbacktrace
-O0
-g
```

## Test results

```text
test_inverse: PASS
test_linalg: PASS
test_portfolio: PASS
test_theta: PASS
validation: PASS
```

The demo and both examples also compile and run in the validation script.

## Reference checks

Fixed reference values were calculated independently with NumPy for a
six-observation, two-asset sample. Tests verify:

- packed unified second moments
- selected entries of empirical covariance
- selected entries of Gaussian covariance
- the complete packed inverse second moment
- the unconstrained Markowitz portfolio
- a one-asset subspace constraint
- a zero-covariance hedging constraint

## Structural checks

Tests also verify:

- `vech`/`ivech` round trips
- duplication-matrix vectorization
- Kronecker-product dimensions and values
- fitted and omitted intercepts
- exact covariance symmetry
- nonnegative covariance diagonals within numerical tolerance
- HAC covariance construction
- conditional coefficient dimensions
- coefficient-covariance extraction
- distinction between upstream and all-column weight modes

## Optimized build

The complete source, tests, demo, and examples are also compiled and run with
`-O2` during release validation.

## FPM

The `fpm` executable was not installed in the validation environment. The
manifest was parsed as TOML and source/application/example/test placement uses
standard FPM automatic target discovery. Direct GNU Fortran validation compiles
the same source graph with stricter diagnostics than a default FPM build.
