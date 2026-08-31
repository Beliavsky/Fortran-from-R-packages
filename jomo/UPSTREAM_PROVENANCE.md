# Upstream provenance

## Source snapshot

Translation source: the supplied `jomo-master` archive corresponding to R
package **jomo 2.7-6** (DESCRIPTION date 2023-04-13, CRAN publication
2023-04-15).

Primary upstream native files reviewed:

- `src/jomo1C.c`
- `src/jomo1ranC.c`
- `src/jomo1ranhrC.c`
- `src/jomo1smcC.c`
- `src/jomo1ransmcC.c`
- `src/jomo1ranhrsmcC.c`
- `src/jomo2comC.c`
- `src/jomo2hrC.c`
- `src/jomo2smcC.c`
- `src/jomo2hrsmcC.c`
- `src/pdflib.c`
- `src/wishart.c`

The R wrappers were used to determine defaults, parameter meanings, model
branches, categorical coding, and the distinction between common and
cluster-specific level-1 covariance models.

## Translation strategy

The native C entry points repeat substantial blocks of MCMC logic for
continuous, categorical, mixed, multilevel, heteroscedastic, and
substantive-model-compatible variants.  The Fortran translation factors those
blocks into reusable typed modules:

- `jomo_single_level`
- `jomo_multilevel`
- `jomo_heteroscedastic`
- `jomo_twolevel`
- `jomo_twolevel_heteroscedastic`
- `jomo_substantive`
- `jomo_smc`, which consolidates the repeated proposal/rebuild algorithms in
  `jomo1smcC.c`, `jomo1ransmcC.c`, `jomo1ranhrsmcC.c`, `jomo2smcC.c`, and
  `jomo2hrsmcC.c`
- shared latent-normal, distribution, RNG, and dense SPD helpers

This is an algorithmic translation, not an ABI translation of R `.Call`
interfaces.  Explicit observation masks replace R `NA`, and one-based
contiguous cluster labels replace factor internals.

## Repository dependency check

Before finalizing the translation, the target `Fortran-from-R-packages`
repository was checked for compatible top-level/shared translations.  It
contains shared `rfortran-core` and `rfortran-linalg` modules and translated
packages including `MASS`, `lme4`, and `survival`; no existing top-level
`jomo` translation was found.  The current `jomo` numerical kernels do not
require those R modeling-package APIs because callers provide design matrices
directly.  The package-local SPD routines are deliberately small, dense, and
self-contained and do not embed BLAS/LAPACK source or system links.
