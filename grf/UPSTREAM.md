# Upstream provenance

- R package: `grf`
- Upstream version: `2.6.1`
- Upstream title: *Generalized Random Forests*
- Upstream URL: https://github.com/grf-labs/grf
- CRAN publication date recorded in `DESCRIPTION`: 2026-03-04
- Upstream license field: `GPL-3`
- Original package tree: `upstream/`

The translation was made from the supplied `grf-master.zip`. The original R and C++ implementation remains under `upstream/`; it is not compiled or vendored into the Fortran library.

## Dependency audit

The target `Fortran-from-R-packages` repository contains translations of `ranger` and `randomForest`, but their public APIs do not supply GRF's honest-tree relabeling and causal/IV split machinery. Reusing either as a dependency would therefore substitute a different estimator rather than reuse a compatible implementation. The maintained Fortran code uses a package-specific tree kernel and small package-local linear algebra routines, and has no external numerical library dependency.

## Translation approach

The Fortran implementation retains the central generalized-random-forest mechanics: randomized feature selection, subsampling, honesty, missing-value routing, family-specific node relabeling, split balance constraints, honest leaf repopulation/pruning, adaptive forest weights, and local moment prediction. It intentionally omits R formula/data-frame/S3 object construction, Rcpp transport, C++ threading, automatic tuning, confidence-interval tree grouping, and presentation functions.

The supplied source package contained the generated file `build/vignette.rds`. It was verified against the upstream MD5 manifest and then omitted from the translation archive to comply with the no-build-products requirement. The original `MD5` manifest itself is retained as provenance.
