# License, attribution, and provenance

This package translates monobin 0.2.4, **Monotonic Binning for Credit Rating
Models**, by **Andrija Djurovic**, author and maintainer (djandrija@gmail.com).
Original CRAN publication: 2022-07-21; upstream license: GPL (>= 3).

- Upstream: https://github.com/andrija-djurovic/monobin
- CRAN: https://cran.r-project.org/package=monobin
- Source: https://cran.r-project.org/src/contrib/monobin_0.2.4.tar.gz
- Source archive SHA-256:
  `8a701a9475fe72ed169efee5a9597571d3e50acdb689ebe8b68cfc26de61c080`
- Retrieved 2026-10-12. Original R sources, DESCRIPTION, NAMESPACE and NEWS are
  retained unmodified under `upstream/`, including their original attribution.
- Fortran adaptation and deterministic test work: Fortran-from-R-packages
  contributors, 2026. Distributed under GPL-3.0-or-later; see `LICENSE`.

The narrow automatic `cut2(g)` adapter follows the algorithm of **Frank E.
Harrell Jr.**, Hmisc `R/cut2.s`: initial source note 3 December 1990, subsequent
modification notes 1992, 1995, 2002 and 2018. Hmisc is authored principally by
Frank E. Harrell Jr., with Cole Beck, Charles Dupont and other contributors, and
licensed GPL (>= 2). GPL-3.0-or-later is used here, compatible with that grant.
No Hmisc R or Fortran dependency source is copied into this package.

Reference automatic grouping was checked against Hmisc 5.3-0:
https://cran.r-project.org/src/contrib/Hmisc_5.3-0.tar.gz
(SHA-256 `11005db1266040f864de3f7f16ac7b533cce2a57a47aad14806d7fff686ca60d`).
The existing top-level Hmisc package is reused for explicit-cut reassignment.

Dependency inspection used Fortran-from-R-packages repository revision
`cbd437d6d92f5481706fd8ea5a7829aaff7b6e24`. Compatible shared APIs are reused:
`r_kinds::dp`, `r_ordering::r_order`, `r_sorting::r_average_ranks`, and
`r_distributions::{r_pnorm,r_pt}` from MIT-licensed rfortran-core, and
`hmisc::cut2_explicit` from GPL-2.0-or-later Hmisc. Hmisc's existing automatic
`cut2(g)` API is absent, requiring only the focused adapter. No shared package
was changed. No BLAS/LAPACK/ARPACK, r.f90, r_mod.f90, or numerical dependency
source is included, and no system-library links are requested.

This translation is not endorsed by the upstream author. README documents
interface changes, preserved quirks, tested comparisons, and numerical limits.
