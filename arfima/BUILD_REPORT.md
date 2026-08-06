# Build report

## Environment

- Compiler: GNU Fortran (Debian 14.2.0-19) 14.2.0
- Fortran standard: Fortran 2018
- FPM: not installed
- Platform: Linux x86-64

## Checked build

Command:

```sh
make check
```

Result: all six test programs passed with warnings treated as errors, runtime
bounds/allocation checks, backtraces, floating-point traps, and uninitialized
real initialization.

## Optimized build

Command:

```sh
make optimized
```

Result: all six test programs passed with `-O3` and warnings treated as
errors.

## Example

Command:

```sh
make example
```

Result: the demonstration simulated, fitted, and forecast an ARFIMA model and
completed successfully with finite estimates and forecast variances.

## FPM validation

The FPM executable was not installed in the validation environment. The
`fpm.toml` file was parsed as TOML and its package, build, install, example,
and test declarations were checked. The same ordered source graph was compiled
and tested through the included Makefile.
