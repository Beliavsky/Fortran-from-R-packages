# Build and translation report

## Source coverage

- 41 fixed-form Fortran source files from the package numerical core were converted to free-form `.f90`.
- The converted core contains approximately 29,600 source lines.
- All 40 R exports are mapped to their underlying numerical entry points in `API_MAP.md`.
- R plotting, S3 methods, and time-series metadata are not part of the standalone Fortran library.

## Modernization applied

- Free-form Fortran source and `.f90` extensions.
- `implicit none` throughout the converted numerical source.
- Shared real kind `dp = kind(1.0d0)` from `timsac_kinds`.
- `real(dp)` declarations in place of `double precision`.
- A module-based public API with allocatable result types for selected common routines.
- A static-library build, an example program, and regression tests.

## Verification

The project was built with GNU Fortran using:

```text
gfortran -std=f2018 -O2 -ffree-line-length-none
```

The regression test executable passed tests for:

- univariate autocorrelation;
- direct matrix filtering;
- recursive matrix filtering; and
- white-noise generation (smoke test).

## Remaining legacy internals

The full numerical core is buildable as free-form Fortran 2018 source, but many original algorithms still use historical constructs such as labeled `do`, `goto`, and `common`. Replacing every internal control-flow and shared-storage construct would be a separate algorithm-by-algorithm refactor and validation project. The current package therefore preserves numerical behavior while providing a modern build and public module layer.
