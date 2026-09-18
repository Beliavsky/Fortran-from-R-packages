# Provenance

## Upstream

- R package: `stinepack`
- Version: 1.5
- Package date: 2024-03-07
- Repository field: CRAN
- Declared license: GPL (>= 2)
- Input archive for this translation: `stinepack-master.zip`

Retained source material is stored under `upstream/`:

- `upstream/DESCRIPTION`
- `upstream/NAMESPACE`
- `upstream/R/na.stinterp.R`
- `upstream/R/parabolaSlopes.R`
- `upstream/R/stinemanSlopes.R`
- `upstream/R/stinterp.R`

The retained R files are included for source traceability and coverage auditing.
They are not compiled by FPM and are not copied dependency source.

## Translation scope

The port translates the numerical computations for Stineman/parabola slope
estimation, rational interpolation, and vector/matrix NaN gap filling. Plotting,
interactive behavior, S3/class metadata, and other R-specific object/interface
facilities are outside the translation scope.

## Dependencies

No external numerical dependency is required. Before implementation, the target
`Beliavsky/Fortran-from-R-packages` repository was checked for an existing
`stinepack` translation or a shared Stineman-interpolation API; none was found.
The implementation therefore uses only standard Fortran intrinsic facilities.

## Working precision

`src/stinepack_kinds.f90` defines the single package real kind `dp = real64`.
All maintained real computations use `real(dp)` and `_dp` constants.
