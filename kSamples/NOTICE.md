# NOTICE and provenance

This repository is a modern free-form Fortran translation of computational code from the R package **kSamples**, version **1.2-12** (2025-08-25), authored by Fritz Scholz and Angie Zhu and distributed by CRAN under **GPL (>= 2)**. The translation was prepared from the user-supplied `kSamples-master.zip`. The original package metadata and the computational R/C sources used for translation are retained under `upstream/` for provenance only and are not built by FPM.

The Fortran translation is distributed under **GPL-2.0-or-later**, consistent with the upstream package. See `LICENSE` for the GNU General Public License version 2 text. Upstream copyright, authorship, attribution, and provenance remain with the original authors and contributors.

## External sibling dependencies

The translation intentionally does **not** vendor dependency source code.

- `rfortran-core` is used for the common `dp` kind, R-like distributions, average ranks, and sample variance helpers.
- `SuppDists` / FPM package `suppdists-fortran` is reused for the normal-order scores corresponding to R `SuppDists::normOrder`.

Both are expected as sibling top-level directories in `Fortran-from-R-packages` and are referenced by FPM path dependencies. Their own licenses and attribution remain in those repositories.

No BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, or dependency implementation is copied into this package. No system `-lblas` or `-llapack` link is required.

## Deliberately omitted R-only surface

`pp.kSamples`, `print.kSamples`, plotting/graphics helpers, formula/list dispatch, presentation formatting, and S3 construction are not translated. These are excluded from the exported computational-function coverage denominator under the translation convention requested for this repository.
