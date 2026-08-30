# NOTICE and licensing

This directory is a modern free-form Fortran translation of computational
parts of the R package **spdep** version **1.4-2** (2026-02-13).

Upstream project:

- Package: `spdep`
- Repository: <https://github.com/r-spatial/spdep/>
- Documentation: <https://r-spatial.github.io/spdep/>
- Upstream maintainer: Roger Bivand
- Upstream license declaration: `GPL (>= 2)`

The upstream author and contributor list is preserved verbatim in
`UPSTREAM_DESCRIPTION`, the upstream citation information is preserved in
`UPSTREAM_CITATION.R`, and file-level upstream copyright lines are retained in
`UPSTREAM_COPYRIGHT_NOTICES.md`. The original package's computational R and C sources
were used as the primary behavioral reference for this translation.

The Fortran translation is distributed under the same **GNU General Public
License, version 2 or (at your option) any later version**. The GPLv2 license
text is in `LICENSE` and `LICENSES/GPL-2.0.txt`. Nothing in this translation is
intended to remove or narrow upstream copyright, attribution, citation, or
license obligations.

No BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, translated R-package dependency,
or other third-party source tree is copied or vendored in this package.
