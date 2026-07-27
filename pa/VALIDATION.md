# Validation

## Toolchain

- GNU Fortran 14.2.0
- LAPACK and BLAS from the system toolchain

## Debug flags

```text
-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror
-ffree-line-length-none -O0 -g -fcheck=all -fbacktrace
```

## Release flags

```text
-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror
-ffree-line-length-none -O2 -fbacktrace
```

## Tested numerical paths

- Average ranks with ties and R-style quintile assignment.
- Mixed design matrices with two categorical variables and one numeric variable.
- Single- and multi-period categorical/continuous exposures.
- Exact single-period Brinson category weights, category returns, quadrant returns, and effects.
- Multi-period Brinson arithmetic, geometric, and linking aggregations.
- Full-rank regression coefficient recovery.
- Active exposures, contributions, and residual attribution.
- Multi-period regression arithmetic, geometric, and linking aggregations.
- CSV Brinson and regression workflows.
- License-header audit.

The tests use hand-derived values and independently calculated multi-period reference values. Additivity identities are checked for category effects and linked aggregate effects.

## Scope statement

All non-plotting computational routines found in the attached package are represented. R S4 classes, formula/data-frame processing, plotting, formatted display methods, and metadata are excluded.
