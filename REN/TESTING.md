# Testing

Preferred:

```text
fpm test
```

Strict GNU Fortran validation without FPM:

```text
./scripts/build_checked.sh
./scripts/build_optimized.sh
```

The checked script compiles the vendored dependencies and REN with Fortran 2018,
runtime bounds checking, backtraces, warnings, and warnings-as-errors, then runs
all tests and smoke-tests both examples.
