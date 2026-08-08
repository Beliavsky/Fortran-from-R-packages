# Validation

Validated with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

The strict validation script compiles the library modules independently, then compiles and runs every test and example.

Test programs:

1. `test_utils` - domain clipping and population sorting.
2. `test_benchmarks` - direct numerical checks of both upstream benchmark functions.
3. `test_minimize` - two-dimensional quadratic minimization.
4. `test_maximize` - one-dimensional maximization and internal score-sign convention.
5. `test_chaos` - bounded elite/bookkeeping behavior with chaos-related options active.

Both examples also compile and execute:

- `quadratic_example`
- `rceim_benchmark_example`

The objective is invoked only through module-level procedures with explicit non-optional procedure dummy interfaces. This design specifically avoids the host-associated optional-callback implicit-interface issue seen with some Windows gfortran/FPM combinations.

FPM itself was not installed in the validation container. `fpm.toml` was syntax-checked with a TOML parser, and the FPM source/test/example layout was independently compiled with the strict scripts.
