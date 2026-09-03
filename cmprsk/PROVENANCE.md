# Provenance

## Upstream

The source archive supplied for this translation identifies:

- Package: `cmprsk`
- Version: 2.2-12
- Date: 2024-05-14
- Author/Maintainer: Robert Gray
- License: GPL (>= 2)
- Native-source copyright: Copyright (C) 2000 Robert Gray

`LICENSE` is copied verbatim from upstream `COPYING` (GNU GPL Version 2). The upstream `DESCRIPTION` is retained under `upstream/DESCRIPTION`.

## Translation mapping

- `src/cincsub.f` -> `src/cmprsk_cuminc.f90` (`cumulative_incidence`).
- `src/crstm.f` -> `src/cmprsk_cuminc.f90` (`gray_test` and its stratum kernel).
- `src/tpoi.f` -> `src/cmprsk_cuminc.f90` (`timepoint_indices`).
- `src/crr.f` -> `src/cmprsk_crr_kernels.f90` (objective, derivatives, variance, residual, and baseline-jump kernels).
- `R/cmprsk.R::cuminc` -> `src/cmprsk_api.f90` (`fit_cuminc` orchestration).
- `R/cmprsk.R::crr` -> `src/cmprsk_crr.f90` (sorting, censoring KM, Newton/Armijo fitting, result assembly).
- `R/cmprsk.R::predict.crr` -> `src/cmprsk_crr.f90` (`predict_crr`).
- `R/cmprsk.R::summary.crr` -> `src/cmprsk_summary.f90` (`summarize_crr` numerical calculations).

The original fixed-form source is not embedded in the distributable package; this keeps all maintained Fortran free-form. Hashes of the supplied upstream computational files are recorded in `UPSTREAM_SHA256SUMS.txt`.

## Shared repository dependencies

Before translation, the target `Fortran-from-R-packages` repository was checked. No existing top-level `cmprsk` translation was present. The translation therefore reuses:

- `rfortran-core` for `dp`, chi-square/normal distribution functions, and shared R-like numerical infrastructure.
- `rfortran-linalg` for checked linear solves and matrix inversion, backed by the repository's pinned pure-Fortran LAPACK dependency.

No source from either dependency is copied into this package.
