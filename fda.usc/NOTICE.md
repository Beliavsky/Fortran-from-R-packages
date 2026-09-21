# NOTICE and provenance

This package is a modern Fortran translation of computational code from the R package **fda.usc 2.2.0** (2024-11-04).

Upstream package authors and contributors, as recorded in `DESCRIPTION`:

- Manuel Febrero Bande - author
- Manuel Oviedo de la Fuente - author and maintainer
- Pedro Galeano - contributor
- Alicia Nieto - contributor
- Eduardo Garcia-Portugues - contributor

Upstream project: `https://github.com/moviedo5/fda.usc`

The upstream package declares the license `GPL-2`. The GNU General Public License version 2 text is included as `COPYING`. This translation preserves that upstream licensing declaration; no relicensing is intended.

## Compiled-code provenance

The upstream files `src/Adot.f90` and `src/PCvM_statistic.f90` state that they were created by **Eduardo Garcia-Portugues**. Their numerical algorithms are translated into `src/fdausc_statistics.f90` as `adot_matrix_vector` and `pcvm_statistic_value`, respectively, using free-form modern Fortran, the package-wide `real(dp)` kind, explicit interfaces, and documented dummy arguments.

The upstream `src/rp_stat.f90` random-projection CvM/KS statistic is likewise translated into `rp_projection_statistics` in `src/fdausc_statistics.f90`.

The R sources remain the authoritative upstream reference for all other translated numerical operations. Exact repository-relative R source paths are recorded in every `[[extra.translation.function]]` entry in `fpm.toml`.

## Shared dependencies

This package does not copy BLAS, LAPACK, or translated dependency sources. It reuses sibling FPM packages:

- `rfortran-core` for `r_kinds::dp` and shared statistical helpers.
- `rfortran-linalg` for the high-level solve, inverse, eigendecomposition, and SVD API.

The target `Fortran-from-R-packages` repository also contains translations of several packages listed by upstream fda.usc, but they are not added as dependencies unless the current Fortran API actually uses them.
