# RobStatTM modern Fortran

This package is a modern Fortran 2018 translation of the numerical core of the R package **RobStatTM 1.0.11**. It is organized as an FPM library and preserves the upstream GPL-3.0-or-later licensing.

The port focuses on reusable numerical procedures rather than R runtime infrastructure. Observations are rows and variables are columns. Regression routines accept an explicit design matrix, so an intercept is included by placing a column of ones in that matrix.

## Implemented numerical areas

- Robust loss functions, derivatives, Gaussian-efficiency tuning, M-scales, and inverse robust R-squared.
- Robust location and scale estimation.
- MM/M regression, DCML mixing and covariance, least-absolute-residual initialization, RFPE, robust nested-model tests, and numeric RFPE stepwise selection.
- BY, leverage-weighted BY, and leverage-weighted maximum-likelihood logistic regression.
- Classical covariance, fast-MVE, projection initialization, MM-SHR scatter, Rocke scatter, Mahalanobis distances, and correlations.
- Residual M-scale robust PCA and a `prcompRob`-compatible numeric interface.

The top-level `robstattm` module exports both idiomatic snake-case procedures and compatibility names such as `lmrobdetmm`, `lmrobdetdcml`, `bylogreg`, `multirobu`, `rockemulti`, `smpca`, and `prcomprob`.

## Build with FPM

The package has local path dependencies under `vendor/` and links LAPACK and BLAS.

```text
fpm test --profile debug
fpm test --profile release
fpm run --example robust_regression_demo
```

## Build without FPM

On a Unix-like system with GNU Fortran, LAPACK, and BLAS:

```text
./scripts/build_checked.sh
./scripts/build_optimized.sh
./build/checked/bin/robust_regression_demo
```

The checked build uses Fortran 2018 conformance, all common warnings as errors, implicit-interface checks, array/bounds/runtime checks, and backtraces.

## Translation boundaries

Plotting, formula parsing, model frames, S3/S4 classes, print/summary methods, R-specific missing-value dispatch, and package datasets are not exposed as Fortran APIs. The complete upstream source snapshot is retained under `upstream/` for provenance.

Several upstream methods depend on `pyinit`, `robustbase`, and `rrcov`. The supplied robustbase Fortran translation is vendored. The previously translated rrcov package is also vendored for MVE, robust PCA initialization, probability functions, and linear algebra. Where the R implementation delegates to Pen~a-Yohai candidate generation or rich formula metadata, this port provides a self-contained numerical equivalent and records the distinction in `TRANSLATION_COVERAGE.md` and `PORTING_NOTES.md`.

## Licensing

The combined package is distributed under GPL-3.0-or-later. The vendored robustbase translation is GPL-2.0-or-later and is used under GPL version 3; the vendored rrcov translation is GPL-3.0-or-later. Full license texts are in `LICENSES/`.
