# NOTICE and provenance

## Upstream project

This package is a modern Fortran translation of computational portions of the R package **psych**, version 2.6.5 (2026-05-11), authored and maintained upstream by William Revelle and other contributors credited by the upstream project. The upstream DESCRIPTION declares `License: GPL (>= 2)`.

Upstream DESCRIPTION, NAMESPACE, and R source files are retained under `upstream/` to preserve attribution and make each translation-coverage mapping auditable. The translated Fortran code is a new implementation and does not copy dependency source code.

## Shared Fortran dependencies

The package uses sibling FPM dependencies `rfortran-core` and `rfortran-linalg` from the Fortran-from-R-packages repository. No source from those packages is embedded here. Linear algebra uses the public high-level `rfortran-linalg` API; this package does not directly link system BLAS/LAPACK and does not directly depend on or vendor fortran-lapack.

The upstream R package imports `mnormt` for some latent-normal probability calculations. This translation instead evaluates the pairwise bivariate-normal probabilities used by its tetrachoric/polychoric routines directly, so no translated `mnormt` dependency is required for the mapped scope.

The upstream `glb.algebraic` function uses `Rcsdp`. This translation does not vendor or depend on Rcsdp; it solves the function's bounded minimum-trace semidefinite problem with a package-local ADMM method and the existing shared eigensolver API. Upstream solver status codes and iteration paths are therefore not reproduced exactly.

## Compatibility notes

R formulas, S3 classes and methods, plotting/graphics, interactive functions, external package orchestration, parallel execution, exact warning text, and R RNG streams are generally outside the typed numerical Fortran API unless explicitly documented otherwise. See README.md and docs/API_COVERAGE.md.
