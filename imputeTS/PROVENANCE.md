# Provenance

## Source snapshot

Translation source: user-supplied `imputeTS-master.zip`, whose `DESCRIPTION` reports package version 3.4 and publication date 2025-08-25. The relevant original files are retained under `upstream/`:

- `DESCRIPTION`, `NAMESPACE`
- computational `R/na_*.R`, `R/statsNA.R`, and `R/internal_algorithm_interface.R`
- computational portion of `R/deprecated_defunct.R`
- native kernels `src/locf.cpp` and `src/ma.cpp`
- `inst/CITATION`

The retained deprecated source file intentionally stops before the upstream plotting-wrapper section. No plotting implementation is translated.

## Dependency review

Before implementing external numerical facilities, the target `Beliavsky/Fortran-from-R-packages` repository was checked. It already contains a top-level `forecast` translation with public ARIMA, frequency-detection, and STL decomposition APIs. Those facilities are referenced through `forecast-fortran = { path = "../forecast" }`.

Stineman interpolation is referenced through `stinepack = { path = "../stinepack" }`, reusing the sibling modern-Fortran translation produced immediately before this package rather than copying its source.

No BLAS/LAPACK system link is introduced by this package and no dependency source is vendored.

## Translation choices

- `src/imputets_basic.f90` translates the direct imputation kernels, including the C++ LOCF and moving-average behavior.
- `src/imputets_interpolation.f90` provides linear and cubic spline interpolation and delegates Stineman interpolation to `stinepack`.
- `src/imputets_kalman.f90` provides a local-linear-trend Kalman smoother and a shared-`forecast` ARMA path. This is the principal intentionally partial mapping.
- `src/imputets_seasonal.f90` delegates period detection and STL to `forecast` and retains the upstream seasonal-decomposition/split orchestration.
- `src/imputets_stats.f90` translates computational missing-gap statistics while omitting printing.

See `API_COVERAGE.md` for function-by-function scope.
