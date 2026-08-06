# Build report

Validated with GNU Fortran 14.2 using:

- Checked mode: Fortran 2018, pedantic diagnostics, warnings as errors, conversion diagnostics, `-fcheck=all`, and backtraces.
- Optimized mode: the same diagnostics with `-O3`.

Five test executables and the demonstration passed in both configurations. The FPM executable was not installed in the validation environment; `fpm.toml` was parsed successfully as TOML and the same source graph was built through the included Makefile.
