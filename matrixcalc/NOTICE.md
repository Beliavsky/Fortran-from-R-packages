# NOTICE and provenance

This project is a modern Fortran translation of the computational code in the R package **matrixcalc** version 1.0-6.

Upstream metadata identifies:

- Author: Frederick Novomestky <fnovomes@poly.edu>
- Maintainer of the referenced release: S. Thomas Kelly <tomkellygenetics@gmail.com>
- Upstream license: GPL (>= 2)
- Upstream bug tracker: https://github.com/TomKellyGenetics/matrixcalc/issues
- CRAN publication date recorded in `DESCRIPTION`: 2022-09-14

The unmodified upstream `DESCRIPTION`, `NAMESPACE`, and R sources used for this translation are retained under `upstream/`.

The Fortran translation is distributed under GPL-2.0-or-later, consistent with the upstream GPL (>= 2) license. The `LICENSE` file contains the GNU General Public License version 2 text.

## Shared Fortran dependencies

The package reuses sibling repository packages rather than vendoring their source:

- `rfortran-core` for the shared `dp` real kind through `r_kinds`.
- `rfortran-linalg` for linear-system solution and symmetric eigendecomposition.

No BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, or translated dependency source is copied into this package. The manifest does not request system `-lblas` or `-llapack` links.
