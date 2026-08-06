# Build report

Validation environment:

- GNU Fortran 14.2.0
- Fortran 2018 mode
- Linux x86-64

Successful commands:

```text
make MODE=checked clean test
make MODE=checked example
make MODE=optimized clean test
make MODE=optimized example
```

All five test programs passed in checked and optimized builds. The checked
build enabled `-fcheck=all`, backtraces, and warning diagnostics. The example
selected a positive penalty and produced a finite periodogram-based standard
error estimate.

FPM was not installed in the validation environment. `fpm.toml` follows the
same source/test/example layout exercised by the Makefile.
