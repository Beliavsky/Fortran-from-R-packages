# Build report

- Compiler: GNU Fortran 14.2.0
- Language mode: Fortran 2018
- Checked build: passed all five tests and the example
- Optimized build: passed all five tests and the example
- External numerical libraries: none
- FPM executable: unavailable in the validation environment; `fpm.toml` was parsed as TOML and the equivalent source graph was built with GNU Make

The checked flags include warnings as errors, bounds/runtime checks,
floating-point traps, backtraces, and signaling-NaN initialization.
