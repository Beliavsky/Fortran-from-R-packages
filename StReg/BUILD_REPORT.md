# Build report

Validation environment:

- GNU Fortran 14.2.0
- Debian Linux
- Fortran 2018 mode

Commands:

```sh
make clean test
make MODE=optimized clean test
make MODE=optimized example
```

All five test programs passed in checked and optimized builds without compiler
or linker warnings. The example converged and recovered an intercept of
approximately 1.20 and a slope of approximately 2.30.
