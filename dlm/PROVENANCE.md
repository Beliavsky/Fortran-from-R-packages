# Provenance

This package translates computational portions of the R package **dlm 1.1-6.1**.
The upstream `DESCRIPTION`, `NAMESPACE`, `R/DLM.R`, `R/arms.R`, `src/dlm.c`, `src/arms-R.c`, and `inst/CITATION`
files are retained verbatim under `provenance/upstream/` as reviewable source
provenance. They are not part of the FPM source tree.

The upstream package declares `License: GPL (>= 2)`. The maintained Fortran
source therefore carries `SPDX-License-Identifier: GPL-2.0-or-later`.

The translation uses the repository's shared `rfortran-core` module for the
single `dp` real kind and `rfortran-linalg` for checked dense linear algebra.
The latter in turn uses the repository-standard pinned pure-Fortran LAPACK FPM
dependency; this package does not request system BLAS/LAPACK libraries. The
`dlmMLE` translation reuses the repository's existing top-level `optimx` package
through a sibling FPM path dependency rather than copying optimizer source.

Numerical algorithms translated here include model construction, stationary AR
parameter transformation, Kalman filtering with missing observations,
likelihood evaluation, Rauch-Tung-Striebel smoothing, deterministic forecasting,
backward state sampling, the `dlmGibbsDIG` precision Gibbs sampler, Wishart
simulation, convex-support boundary search, the Gilks adaptive rejection
Metropolis sampler and multivariate hit-and-run wrapper, `dlmMLE` likelihood
optimization orchestration, and MCMC output diagnostics.
See `API_COVERAGE.md` and the `[extra.translation]` manifest tables for exact
function-level claims and omissions.

## Source archive integrity

The uploaded upstream archive used for this translation has SHA-256:

`22b4edfb0d79c65b869d7ca1ea5d0f964b713aa643459ede6c5213c94e89b4d5`

Per-file hashes for the retained review sources are recorded in
`provenance/SHA256SUMS`.
