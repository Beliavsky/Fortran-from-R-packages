# Validation

Compiler: GNU Fortran 14.2.0.

## Strict build

All nine regression executables pass with:

```text
-std=f2018 -O0 -g -fcheck=all
-Wall -Wextra -Werror
-Wimplicit-interface -Werror=implicit-interface
```

Tests:

1. `test_mechanisms` -- Boston/DA, blocking pairs, RSD, tenant TTC.
2. `test_eadam` -- no-consent EADAM equality with ordinary deferred acceptance.
3. `test_ttc_variants` -- school TTC and TTC-and-chains.
4. `test_stable_sets` -- exact HRI and SRI stable-solution enumeration.
5. `test_couples` -- joint couples assignment and stability fallback.
6. `test_plp` -- exact binary partitioning LP and objective.
7. `test_stats` -- matrix-level KHB/probit/OLS calculations.
8. `test_sim` -- one- and two-sided simulation dimensions/capacities.
9. `test_helpers` -- pair/coalition construction and consensus Monte Carlo.

Both examples compile and run in the same build.

## Optimized build

The same nine tests and both examples pass using `-O2` with all warnings and
implicit-interface diagnostics promoted to errors.

## FPM metadata

FPM itself was not installed in the validation container. The top-level and
both vendored dependency `fpm.toml` files were parsed successfully with a TOML
parser, and the source layout was compiled using the equivalent direct
`gfortran` dependency order.

## Portability checks

- No main-package Fortran source/test/example line exceeds 132 columns.
- Every external LAPACK procedure has an explicit interface.
- The code does not rely on short-circuit `.and.` / `.or.` evaluation for
  array-index safety.
