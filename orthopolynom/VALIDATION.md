# Validation

The port was validated with GNU Fortran in two configurations.

## Runtime-checked build

Compiler options:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -fbacktrace
```

Results:

```text
test_families:   PASS
test_helpers:    PASS
test_parametric: PASS
```

## Optimized build

Compiler options:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface
```

Results:

```text
test_families:   PASS
test_helpers:    PASS
test_parametric: PASS
```

The example program also runs and returns the expected degree-five Legendre
coefficients and Chebyshev-T roots.

Independent parameterized-family reference values were generated from standard
special-function definitions for Gegenbauer, Jacobi P, generalized Laguerre,
monic shifted Jacobi (the package's Jacobi G convention), and generalized
Hermite polynomials.

All 88 names exported by the upstream NAMESPACE have a mapped computational
Fortran implementation. Source, test, and example files stay within the
standard 132-column free-form Fortran line limit.

FPM itself was not installed in the validation environment. `fpm.toml` was
parsed as TOML and the exact `src/`, `test/`, `example/`, and local dependency
layout was compiled directly with GNU Fortran.
