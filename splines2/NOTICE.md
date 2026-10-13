# Provenance, copyright, licenses, dependencies

## Upstream package

- R package: **splines2 0.5.4**, CRAN publication 2025-02-28.
- Authors: **Wenjie Wang and Jun Yan**.
- Copyright **(C) 2016-2025**, as retained in upstream R sources and numerical
  C++ headers; translated maintained source retains this notice.
- License: **GNU GPL version 3 or any later version**. The complete GPLv3 text
  is in `LICENSE`; this package is distributed under GPL-3.0-or-later.
- Pinned CRAN mirror: https://github.com/cran/splines2/tree/2321fc1788c387ad06e531d408c31b5d840b2286
- Original project: https://github.com/wenjie2wang/splines2
- Package site: https://wwenjie.org/splines2
- Algorithm reference: Wang and Yan (2021), https://doi.org/10.6339/21-JDS1020

The original `DESCRIPTION`, `NAMESPACE`, and applicable exported R source files
are retained under `upstream/`, including their full file notices. The numerical
translation was checked against `src/splines2_export.cpp` and these original
headers under `inst/include/splines2Armadillo/`: `SplineBase.h`, `BSpline.h`,
`MSpline.h`, `ISpline.h`, `CSpline.h`, `NaturalSpline.h`, `NaturalSplineK.h`,
`PeriodicSpline.h`, `BernsteinPoly.h`, and `utils.h`. Their algorithms include
normalizations, transforms, knot extension and derivative/integral recurrences.
No Rcpp or Armadillo source/library is copied or required by the translation.

## Existing repository dependencies inspected and reused

Repository inspected at commit **cbd437d6d92f5481706fd8ea5a7829aaff7b6e24** of
https://github.com/Beliavsky/Fortran-from-R-packages .

- `splines` / FPM name `splines-fortran`: GPL-2.0-or-later, copyright (C) 1998
  Douglas M. Bates and William N. Venables; modern Fortran translation 2026.
  Reuses public `splines_core::spline_basis_nonzero` and
  `splines_linalg::solve_linear`. GPL-2.0-or-later permits GPLv3-compatible use.
  The existing solver is sufficient for the small cardinal-basis transformation,
  and does not require a new LAPACK dependency or duplicate solver.
- `rfortran-core`: MIT, Vivek Rao and contributors. Reuses `r_kinds::dp`,
  `r_quantiles::r_quantile_type7`, and `r_ordering::r_sort_values_in_place`.
- `pbs` was inspected, including its pure B-spline evaluator. It declares
  GPL-2.0-only, so it is not used by this GPL-3.0-or-later package.
- `rfortran-linalg` was inspected. No new BLAS/LAPACK layer is needed because
  the existing compatible spline package already provides the required solve.

Dependencies remain sibling FPM paths. Their source files, licenses, binaries,
R compatibility modules, BLAS, LAPACK, and other translated packages are not
vendored here. No system BLAS/LAPACK links or runtime R dependencies occur in the
manifest. Dependency notices remain with their separately maintained packages.

## Translation and tests

Modern free-form Fortran adaptation and tests, 2026. Maintained package sources
are GPL-3.0-or-later. The deterministic reference dataset is generated from the
pinned upstream package using `tools/generate_reference.R`; the script documents
the deliberate upstream C-spline derivative bug repair. See `VALIDATION.md`.
