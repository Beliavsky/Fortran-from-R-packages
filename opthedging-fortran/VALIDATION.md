# Validation report

Validation date: 2026-07-25

## Compiler

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
```

## Strict debug build

The project was compiled directly because FPM was not installed in the
validation runtime. `scripts/validate_gfortran.sh` uses:

```text
-std=f2018 -O0 -g -Wall -Wextra -Wconversion-extra
-Wimplicit-interface -Werror -fcheck=all -fbacktrace
-ffree-line-length-132
```

Results:

```text
test_hedging: PASS
test_interpolation: PASS
test_parity: PASS
test_reference: PASS
validation: PASS
```

The application and both examples also compiled and ran successfully.

## Optimized build

The same suite was rebuilt with `-O2`; all tests, the application, and examples
passed.

## Test coverage

- Interior interpolation.
- Left and right linear extrapolation.
- One-point interpolation safety.
- Direct one-period comparison with independently evaluated formulas.
- Fixed multi-period reference values computed independently in NumPy.
- Direct comparison with a standalone build of the original C recursion.
- Multi-period discounted put-call parity across the complete grid.
- Call-minus-put hedge parity of one share.
- Grid and interpolated result methods.
- Rejection of a singular constant-return sample.
- Demonstration and example execution.

## Release audits

The release was checked for:

- Valid TOML syntax.
- ASCII-only authored source and documentation.
- Maximum free-form Fortran line length of 132 columns.
- `implicit none` in every Fortran program unit file.
- GPL SPDX identifiers in every Fortran file.
- Original and translated SHA-256 manifests.
- Absence of build objects and module files from the release archive.

A deterministic five-period fixture was also run through a standalone build of
the original C recursion. `rho`, all five option values, and all five auxiliary
coefficients at the central grid point agreed to rounding error at approximately
machine precision.

The exact final ZIP was extracted into a clean directory and rebuilt with the
strict validation script.
