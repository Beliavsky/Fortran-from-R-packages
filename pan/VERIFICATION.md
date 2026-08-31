# Verification

This translation was verified in the available Linux execution environment with GNU Fortran 14.2.0.

## Strict compiler verification

The library modules, all tests, and the example were compiled with:

```text
-std=f2018 -Wall -Wextra -Werror=implicit-interface -fcheck=all
```

Results:

- library compilation: PASS
- deterministic tests: 9/9 PASS
- examples: 1/1 PASS

The full compiler/test transcript is in `GFORTRAN_VERIFICATION.txt`.

## Source-policy audit

Seventeen maintained Fortran source/test/example files were checked. The audit found zero violations for:

- free-form lines longer than 132 characters
- semicolon-separated executable statements
- `double precision`, `real*8`, `kind(0.0d0)`, or D-exponent real literals
- self-comparison NaN tests
- dummy arguments lacking `INTENT`/`VALUE`
- dummy declarations combining multiple dummy arguments
- dummy declarations lacking meaningful trailing FORD `!!` documentation
- duplicate Fortran source files
- copied `r.f90`/`r_mod.f90`
- compiled/build/archive artifacts in the package tree
- system BLAS/LAPACK links or unsafe fast-math flags

## FPM availability limitation

The requested FPM commands were explicitly attempted:

```text
fpm build
fpm test
fpm clean --all
```

The execution environment does not provide an `fpm` executable, so each command returned exit status 127 (`fpm: command not found`). The exact transcript is in `FPM_ATTEMPTS.txt`.

Because FPM was unavailable, successful FPM execution cannot be claimed. The package follows the standard FPM `src/`, `test/`, and `example/` layout and was instead compiled and exercised directly with gfortran using the strict flags above.
