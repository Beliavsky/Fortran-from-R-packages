# Validation

Validation was performed on 2026-09-15 with GNU Fortran 14.2.0 using the maintained source files named by `fpm.toml`.

## Strict debug/runtime-check configuration

Compiler flags:

```text
-std=f2018 -O0 -fcheck=all -fbacktrace -Wall -Wextra -Werror -pedantic
```

Results:

- package sources compile successfully;
- all deterministic tests pass;
- example compiles and runs;
- independent SciPy B-spline parity passes;
- static source/style/coverage audit passes.

Independent parity differences:

```text
periodic max abs diff: 1.110e-16
ordinary max abs diff: 3.997e-15
```

The independent implementation uses SciPy `BSpline` basis functions and derivatives, constructs the periodic wrapped knot vector independently, folds the endpoint basis columns independently, and evaluates the ordinary extrapolation through boundary derivative Taylor terms.

## Optimized configuration

Compiler flags:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -pedantic
```

The full deterministic test suite, example, source audit, and independent parity checks also pass in this configuration.

## FPM and fprettify

FPM is not installed on this execution host. The required commands were explicitly attempted:

```text
fpm build       -> exit 127, command not found
fpm test        -> exit 127, command not found
fpm clean --all -> exit 127, command not found
```

`fprettify` is likewise unavailable (`exit 127`). No successful FPM/fprettify run is claimed. Maintained Fortran source is nevertheless free form, warning-clean under the compiler configurations above, and within the default 132-column free-form line limit.

`VALIDATION_RESULTS.txt` records the command results. Before archive creation the tree is manually checked for compiler products, cache directories, nested ZIP files, duplicate Fortran source, and copied dependency source.
