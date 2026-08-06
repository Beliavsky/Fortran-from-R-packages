# Build report

## Environment

- Compiler: GNU Fortran 14.2.0
- Language mode: Fortran 2018
- Platform used for validation: 64-bit Linux
- External numerical libraries: none
- FPM executable: not available in the validation environment

## Checked build

The checked build uses:

```text
-std=f2018 -Wall -Wextra -Wpedantic -Werror -Wno-compare-reals
-fimplicit-none -ffree-line-length-none -O0 -g -fcheck=all
-fbacktrace -finit-real=snan
```

All five test programs passed, followed by a successful example run.

## Optimized build

The optimized build uses the same language and warning flags with `-O3`.
All five test programs passed without warnings, followed by a successful
example run.

## FPM manifest

`fpm.toml` was parsed with Python's TOML parser and its package metadata was
validated. The equivalent source graph was compiled and tested through the
included Makefile.

## Archive validation

The release procedure removes all build products, creates a source-only ZIP,
extracts it into an independent directory, and reruns checked and optimized
tests plus the example from the extracted copy.
