# Build report

Validated with GNU Fortran 14.2 on a 64-bit Linux environment.

## Checked build

Command:

```text
make clean
make check
```

Result: all six tests passed with Fortran 2018 conformance, warnings as errors,
bounds/runtime checking, and backtraces enabled.

## Optimized build

Command:

```text
make optimized
```

Result: all six tests passed with `-O3 -march=native` and warnings as errors.

## Demonstration

Command:

```text
make demo
```

Result: the example completed with a finite fitted log likelihood, valid fitted
GJR-GARCH coefficients, and positive five-step forecast standard deviations.

## FPM

The package contains an explicit `fpm.toml` for the library, six test programs,
and one example. An FPM executable was not installed in the validation
environment. The manifest was parsed as TOML, and the same source dependency
order was compiled and tested through the included Makefile and scripts.
