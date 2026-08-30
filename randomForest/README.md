# randomForest - modern Fortran computational translation

This directory is a modern free-form Fortran translation of the computational core of the R package **randomForest 4.7-1.2**. It is intended to live as the top-level `randomForest/` directory in `Beliavsky/Fortran-from-R-packages`.

The translation preserves the package's GPL-2-or-later licensing and upstream algorithm provenance while replacing the R/C/legacy-Fortran interfaces with typed Fortran modules and FPM packaging. Plotting, formula/data-frame/S3 presentation, printing, and other R-only interfaces are intentionally omitted.

## Dependencies

The package reuses sibling shared packages rather than copying numerical infrastructure:

```toml
rfortran-core   = { path = "../rfortran-core" }
rfortran-linalg = { path = "../rfortran-linalg" }
```

`dp` is imported from `r_kinds` in `rfortran-core` and re-exported by the public `randomforest` module. Every maintained real declaration uses `real(dp)` and real constants are kinded with `_dp` where appropriate. `rfortran-linalg` is used only for the symmetric eigensystem in classical MDS; its `real64` API is the same real kind represented by `r_kinds::dp`.

## Build and test

From this directory, with `../rfortran-core` and `../rfortran-linalg` present:

```text
fpm build
fpm test
fpm run --example classification_example
fpm run --example regression_example
```

## Main API

```fortran
use randomforest, only : dp, rf_options
use randomforest, only : rf_classification_forest, rf_regression_forest
use randomforest, only : fit_classification, predict_classification
use randomforest, only : fit_regression, predict_regression
```

The public module additionally exports unsupervised forests, proximity imputation, OOB-based `mtry` tuning, random-forest cross-validation feature selection, classical MDS coordinates, margins, outlier scores, class centers, tree-size/variable-use utilities, and numerical partial dependence.

Categorical predictors are passed as integer-valued columns stored in `real(dp)` arrays, together with `ncat(j) > 1`; numeric predictors use `ncat(j) = 1`. The implementation supports up to 53 categorical levels.

See `API_COVERAGE.md` for routine-by-routine scope and `NOTICE.md` / `UPSTREAM.md` for provenance.
