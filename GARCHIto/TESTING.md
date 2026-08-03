# Testing

## FPM

```text
fpm test
fpm run
fpm run --example unified_fit
fpm run --example realized_fit
fpm run --example option_fit
```

## Direct GNU Fortran validation

Linux/macOS shell:

```text
./scripts/test_gfortran.sh
```

Windows command prompt:

```text
scripts\test_gfortran.bat
```

The scripts compile with Fortran 2018 conformance, all common warnings treated
as errors, bounds checking, runtime checks, and floating-point traps.

## Test coverage

- `test_unified`: synthetic unified GARCH-Ito fit and broad parameter recovery
- `test_realized`: no-jump and jump realized models plus invalid negative data
- `test_option`: all four option-data variants
- `test_validation`: length errors and safe return at a one-iteration limit

The repository is also validated with an optimized `-O3 -Werror` build and by
extracting and rebuilding the final source-only ZIP.
