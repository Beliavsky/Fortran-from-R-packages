# Build report

Translation version: 0.1.0

Compiler used for validation:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
```

## Checked build

Flags include Fortran 2018 conformance, warnings, conversion warnings, full
runtime checking, and backtraces. All five tests and the example passed without
compiler warnings.

## Optimized build

Flags include `-O3 -march=native -Wall -Wextra -Wpedantic`. All five tests and
the example passed without compiler warnings.

## FPM

No `fpm` executable was installed in the validation environment. `fpm.toml` was
parsed successfully as TOML and its local dependency paths were verified. The
same module graph was built with the included Makefile.

## Validated example output

```text
Mean estimate:    0.0039851  iid SE:    0.0019398
Prewhitened correlated SE:    0.0131417  AR(1):   0.79812
95% ES estimate:    0.0262688  iid SE:    0.0006520
```
