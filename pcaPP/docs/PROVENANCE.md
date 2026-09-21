# Provenance

## Upstream

- R package: `pcaPP`
- Upstream version in the supplied archive: `2.0-5`
- Upstream authors: Peter Filzmoser, Heinrich Fritz, Klaudius Kalcher, and
  Valentin Todorov
- Upstream license: GPL (>= 3)
- Upstream project URL recorded in DESCRIPTION:
  `https://github.com/valentint/pcaPP`

The supplied archive's `DESCRIPTION`, `NAMESPACE`, `ChangeLog`, README, and R
source files used for function-level coverage are retained under `upstream/`.
They establish provenance and are not compiled by FPM.

## Translation strategy

The translation preserves numerical operations while replacing dynamic R
interfaces with explicit Fortran arguments and result types.  In particular:

- robust centers are represented by `median_result`;
- scaling is represented by `scale_result` and integer method selectors;
- PCA and sparse PCA return `pca_result`;
- covariance reconstruction returns `covariance_result`;
- sparse-PCA tuning scans return `tuning_result`.

The robust PCA algorithms are independent Fortran implementations.  No
upstream C/C++ source is compiled or copied into `src/`.

## Dependency review

The target `Fortran-from-R-packages` repository was checked before choosing
this implementation strategy.  A top-level `mvtnorm` translation exists, but
pcaPP's translated computational surface only needs standard normal draws for
the simulation/random-candidate paths.  Using Fortran's intrinsic RNG plus a
Box-Muller transform avoids an unnecessary dependency.  No system BLAS or
LAPACK link is needed.

## Important compatibility differences

- Random-number streams are Fortran RNG streams, not R RNG streams.
- R function-valued centering/scaling callbacks are replaced by method codes
  or explicit center/scale vectors.
- Optimizer return/status details do not duplicate R `optim`/`nlm` objects.
- `PCAgrid`, `PCAproj`, and `sPCAgrid` return typed numerical results rather
  than R `princomp`/S3 objects.
- The translation prioritizes the core robust projection algorithms.  Some
  upstream high-dimensional acceleration, tracing, and control options are
  represented more simply.
- `cor.fk` is O(n^2) rather than the upstream O(n log n) Kendall algorithm.
- `qn` uses direct pairwise distances rather than the upstream specialized
  selection implementation.
