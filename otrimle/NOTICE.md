# NOTICE and provenance

## Upstream package

- R package: **otrimle**
- Version: **2.0**
- Title: *Robust Model-Based Clustering*
- Authors: Pietro Coretto and Christian Hennig
- Upstream publication date: 2021-05-29
- Upstream license: GPL (>= 2)
- CRAN package sources supplied by the user are retained under `upstream/` for provenance and coverage auditing.
- The upstream `DESCRIPTION`, `NAMESPACE`, R sources, and `inst/CITATION` are preserved without substantive edits.

The package documentation cites:

- Coretto, P. and Hennig, C. (2016), robust improper maximum-likelihood clustering.
- Coretto, P. and Hennig, C. (2017), OTRIMLE tuning/model-selection methodology.

See the preserved upstream citation metadata in `upstream/inst/CITATION` for the exact references.

## Translation provenance

This is an independent modern free-form Fortran translation of the computational portions of `otrimle`.
It is intended for placement as the top-level directory `otrimle` in
`https://github.com/Beliavsky/Fortran-from-R-packages`.

The translation preserves the upstream robust-improper-mixture ECM structure, global eigenratio restriction,
noise-proportion constraint, weighted chi-square discrepancy criterion, log-improper-density search,
density diagnostics, simulation generator, cluster-count grid evaluation, and simulation-based selection workflow.

## Dependency audit

The upstream R package imports `mclust`, `robustbase`, `mvtnorm`, and parallel/foreach infrastructure.

- `InitClust()` now reuses the repository's existing top-level `mclust` translation through
  `mclust-fortran = { path = "../mclust" }`.  No `mclust` source is copied into this package.  The call uses
  `hc_fit` and `hclass`, corresponding to upstream `mclust::hc(..., minclus=G, modelName=...)` and
  `mclust::hclass`.  The upstream Manhattan-average-link fallback is retained locally; only the final randomized
  retry is deterministic because exact RNG parity is outside scope.
- The repository's translated `robustbase` package was checked before implementing the simulation-summary helper.
  Its current `scale_tau2(x, center)` interface returns only a scale and does not provide the `mu.too=TRUE`
  location/scale pair required by `summary.otrimlesimgdens`; its numerical definition also differs from the
  current upstream default.  Therefore `otrimle` contains a focused translation of the required default
  `robustbase::scaleTau2` formula rather than introducing an incompatible dependency.
- Multivariate-normal simulation remains package-local because only the small simulation surface is needed and RNG
  streams are intentionally not required to match R.
- Parallel `foreach`, `mclapply`, and `doParallel` execution is intentionally replaced by serial loops.

## Additional numerical provenance

`kerndensmeasure()` mirrors the Gaussian path of the R `stats::density.default` implementation used in the era of
this 2021 `otrimle` release: `bw.nrd0`, weighted linear binning (`BinDist`), zero-padded Gaussian convolution, and
linear interpolation to the requested grid.  R's density implementation is Copyright the R Core Team under GPL-2
or later; the weighted binning source also records changes by Adrian Baddeley.  This package is GPL-2-or-later and
contains an independent Fortran implementation of that numerical procedure.

The focused `scaleTau2` translation follows `robustbase/R/OGK.R`, whose source credits Kjell Konis and subsequent
robustbase maintenance by Martin Maechler and contributors.  `robustbase` is GPL-2-or-later.

No source from BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, `rfortran-compat`, `mclust`, `robustbase`, or `mvtnorm`
is copied into this package.
