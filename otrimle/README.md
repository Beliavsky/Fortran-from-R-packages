# otrimle — modern Fortran translation

This directory translates the computational core of the R package `otrimle` 2.0 to modern free-form Fortran.
It implements robust improper maximum-likelihood estimation (RIMLE), OTRIMLE tuning over the constant improper
density, the package's eigenratio and noise-proportion constraints, density diagnostics, and the simulation/model-
selection utilities.

The maintained `otrimle` sources use one real kind (`dp = real64`) throughout and do not vendor BLAS, LAPACK,
ARPACK, `rfortran-compat`, or translated-package source.  `InitClust` reuses the repository's sibling `mclust`
translation through an FPM path dependency, matching the upstream package's first-choice hierarchical initializer.

## Build

From this directory in the target repository:

```text
fpm build
fpm test
```

Run the example with:

```text
fpm run --example otrimle_example
```

The source is standard free-form Fortran and is also validated directly with gfortran; see `VALIDATION.md`.

## Public Fortran API

Import the facade module:

```fortran
use otrimle_mod
```

The main fit routines are:

- `init_clust` — kNN trimming, `mclust` model-based hierarchical initialization, and the upstream Manhattan fallback.
- `rimle` — robust improper maximum-likelihood mixture fitting at a fixed/imputed `logicd`.
- `otrimle_fit_grid` — OTRIMLE search over candidate log improper densities.
- `kerndensmeasure`, `kerndensp`, `kerndenscluster` — density-shape diagnostics.
- `generator_otrimle` — simulation from a fitted improper Gaussian mixture.
- `otrimleg` — evaluate candidate numbers of clusters.
- `otrimlesimg` and `summarize_otrimlesimgdens` — simulation-based density selection.
- `kmeanfun`, `ksdfun` — upstream calibration formulas.

Fitted models use the public derived type `otrimle_fit`, containing the improper log likelihood, criterion, mixture
probabilities, means, covariance matrices, posteriors, squared Mahalanobis distances, classification, and component
sizes.

## Translation coverage

Package status: **substantial**.

**12 of 12 (100.0%)** exported computational functions/methods are mapped. The percentage is a function-mapping
measure; it does **not** mean the full R package interface or every numerical implementation detail is identical.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `InitClust` | `otrimle_mod` | `init_clust` | substantial |
| `otrimle` | `otrimle_mod` | `otrimle_fit_grid` | substantial |
| `rimle` | `otrimle_mod` | `rimle` | substantial |
| `kerndensmeasure` | `otrimle_mod` | `kerndensmeasure` | substantial |
| `kmeanfun` | `otrimle_mod` | `kmeanfun` | complete |
| `ksdfun` | `otrimle_mod` | `ksdfun` | complete |
| `kerndensp` | `otrimle_mod` | `kerndensp` | substantial |
| `kerndenscluster` | `otrimle_mod` | `kerndenscluster` | substantial |
| `generator.otrimle` | `otrimle_mod` | `generator_otrimle` | substantial |
| `otrimleg` | `otrimle_mod` | `otrimleg` | substantial |
| `otrimlesimg` | `otrimle_mod` | `otrimlesimg` | substantial |
| `summary.otrimlesimgdens` | `otrimle_mod` | `summarize_otrimlesimgdens` | substantial |

### Compatibility differences

- R data frames, S3 classes, names/dimnames, warnings, printing, and plotting are replaced by explicit numeric arrays
  and Fortran derived types.
- `InitClust()` now follows the upstream hierarchy: kNN trimming, `mclust::hc`/`hclass` via the sibling Fortran
  `mclust` package, then Manhattan average-linkage if that fit is invalid.  The last-resort ten-try R fallback is
  deterministic rather than randomized, because exact R RNG parity is intentionally out of scope.
- Parallel R execution is replaced by serial deterministic loops.
- `kerndensmeasure()` reproduces the R 4.1-era Gaussian `stats::density()` algorithm used by this 2021 upstream
  snapshot: `bw.nrd0`, weighted linear binning, zero-padded convolution, and linear interpolation.  It does not
  reproduce R's NA handling or density-object metadata.
- `summary.otrimlesimgdens()` now uses a focused translation of the upstream default
  `robustbase::scaleTau2(..., mu.too=TRUE)` calculation rather than median/MAD.  The existing repository
  `robustbase` helper was audited but is not used because its current `scale_tau2` API returns only a scalar scale
  and is therefore not compatible with this upstream call.
- Simulation uses the intrinsic Fortran RNG, not R's RNG stream.

These remaining differences are why the package is marked **substantial**, not complete, despite 12/12 mapped
functions.

## Validation

The deterministic Fortran tests exercise the empirical calibration formulas, initialization, fixed-`logicd` RIMLE,
OTRIMLE grid selection, density diagnostics, simulation, cluster-count evaluation, and simulation-summary workflow.
An independent NumPy/SciPy parity script recomputes the fitted mixture density, posterior probabilities, weighted
chi-square discrepancy criterion, covariance eigenratio constraint, R-compatible binned/interpolated density
diagnostic, and `scaleTau2` robust location/scale from exported values.

See `VALIDATION.md` for the commands and current results.

## License and provenance

The upstream package is GPL (>= 2). See `LICENSE.md`, `COPYING`, `NOTICE.md`, and `upstream/`.
