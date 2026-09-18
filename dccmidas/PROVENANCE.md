# Provenance

## Upstream

- Package: `dccmidas`
- Upstream version inspected: `0.1.3`
- Upstream author/maintainer: Vincenzo Candila
- License declared by upstream: GPL-3
- Input archive: `dccmidas-master.zip`

The following upstream files are retained as reference material:

- `upstream/DESCRIPTION`
- `upstream/NAMESPACE`
- `upstream/R/functions.R`
- `upstream/R/RcppExports.R`
- `upstream/src/det.cpp`
- `upstream/src/invf.cpp`
- `upstream/inst/CITATION`
- `upstream/inst/REFERENCES.bib`

## Translation approach

The R time-series/list/S3 interface is separated from the numerical kernels. Return data are supplied to Fortran as dense numeric matrices, standardized residuals as `(assets,time)` arrays, and conditional covariance paths as `(assets,assets,time)` arrays.

The following model families are translated directly from the upstream formulas:

- DCC-MIDAS and asymmetric DCC-MIDAS
- corrected DCC, asymmetric DCC, and DECO
- scalar and diagonal BEKK
- RiskMetrics and moving covariance
- covariance-loss evaluation and QML standard errors

The high-level `dcc_fit` workflow reuses sibling translations instead of embedding their algorithms:

- `rugarch` for sGARCH, gjrGARCH, eGARCH, iGARCH, and csGARCH first-stage fits
- `rumidas` for GM/DAGM GARCH-MIDAS first-stage fits and MIDAS lag weights
- `maxLik` for multivariate second-stage optimization
- `roll` for the rolling weighted long-run correlation calculation
- `rfortran-linalg` for determinant, inverse, and symmetric eigensystem operations

`bekk_fit` likewise reuses `maxLik` rather than embedding an optimizer.

## Preserved upstream conventions

Several formulas that might otherwise look unusual are intentional upstream behavior and are preserved:

- `riskmetrics_mat` weights the lagged return outer product by `lambda` and the previous covariance by `1-lambda`.
- `moving_cov` averages raw lagged return outer products without demeaning each rolling window.
- corrected DCC uses lagged diagonal rescaling before its innovation outer product.
- asymmetric DCC uses the upstream negative-residual covariance intercept correction and its original start-index conventions.
- asymmetric DCC-MIDAS likelihood starts at the upstream time index, while its matrix estimator fills the initial covariance matrices with the later-path mean.
- `cov_eval` preserves the package's exact SFROB, QLIKE, and RMSE formulas rather than substituting similarly named textbook losses.

## Material compatibility differences

- R lists, xts/zoo date alignment, names, S3 objects, print/summary/plot methods, and presentation code are not translated.
- Missing-data merging performed by xts is not reproduced; callers should provide already aligned complete numeric arrays.
- Out-of-sample rolling forecast orchestration in `dcc_fit` and `bekk_fit` is not translated.
- The Fortran `dcc_fit` first stage supports `norm` and `std`; other rugarch distribution strings from the R interface are not currently mapped.
- GARCH-MIDAS first-stage fitting accepts precomputed lag matrices `(K+1,time,assets)` instead of reproducing the R-side date/frequency conversion of the `MV` list.
- Upstream `bekk_fit` selects among `R` random Uniform starting vectors. The Fortran wrapper uses a deterministic default start unless the caller supplies a starting vector.
- Optimizer iteration paths and floating-point rounding can differ between R/maxLik and the translated Fortran `maxLik`, even when they optimize the same translated likelihood.

These differences are why package status is `substantial` rather than `complete`, even though every exported computational R function in the stated coverage basis has a meaningful Fortran mapping.
