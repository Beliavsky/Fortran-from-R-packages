# Notices and provenance

## Upstream matrixTests

- Package: matrixTests 0.2.3.1, published on CRAN 2025-10-17.
- Author and maintainer: Karolis Koncevičius.
- Source: https://cran.r-project.org/src/contrib/matrixTests_0.2.3.1.tar.gz
- Project: https://github.com/karoliskoncevicius/matrixTests
- Verified source tarball SHA-256:
  `ccbf024557331a2ef0156f1c7c9c61f1c07c711c52bb8302aa5dd0d0411c2666`
- Upstream license: GPL-2, meaning GNU General Public License version 2 only.
  The complete license is in `LICENSE`. Upstream DESCRIPTION and NAMESPACE,
  and the unmodified R computational source files, are retained under
  `upstream/` for attribution, inventory, and source-to-source review.
- Modern free-form Fortran translation: Fortran-from-R-packages contributors,
  2026, distributed under GPL-2.0-only. This work is not endorsed by upstream.

## Shared dependency inspection and reuse

No source from these dependency packages is embedded in this package.

- `rfortran-core`, by Vivek Rao and contributors, MIT: reused named `dp` kind,
  descriptive summaries, ranking/ordering, normal/Student/F/chi-square
  distributions and quantiles. These are compatible common numerical APIs.
- `matrixStats`, by Henrik Bengtsson and contributors, Artistic-2.0 and
  GPL-2.0-or-later in the sibling translation: reused row medians for
  Brown-Forsythe/Fligner. GPL-compatible dependency notices remain in that
  sibling; no matrixStats sources are copied here.
- `exactRankTests`, by Torsten Hothorn and Kurt Hornik, GPL-2.0-or-later:
  reused pure exact permutation distributions and inclusive one-sided tails
  for untied Wilcoxon statistics. Its broader tied exact-inference API is not
  substituted for matrixTests's required normal approximation with ties.
- `nortest`, by Juergen Gross and Uwe Ligges, GPL-2.0-or-later: reused its
  compatible pure Anderson-Darling fitted-normal test including minimum n=8
  and its probability floor.
- `moments`, by Lukasz Komsta and Frederick Novomestky, GPL-2.0-or-later:
  reused compatible Jarque-Bera, skewness, and kurtosis APIs. These procedures
  are impure in the existing sibling, so the calling adapters cannot be pure.
- `stats` was inspected: its available KS helper is asymptotic two-sided only,
  and the sibling's current FPM build lacks `r_matrix_statistics`. It cannot
  supply the required exact tied two-sample inference. The local KS kernel
  is an original direct-upper-tail random-label recursion and converged
  limiting-series implementation, not a copy of a dependency.
- `rfortran-linalg` was inspected and its high-level pivoted least-squares API
  smoke-tested. That API is impure. The user's requested pure numerical call
  chain instead uses direct `la_dgelsy` from the FPM-compatible Beliavsky fork.
  Replacing rfortran-linalg's existing transitive LAPACK source with this fork
  is not API-compatible (`la_lapack_d` is absent), so no such override or shared
  sibling modification is used.

## fortran-lapack

Original modern fortran-lapack author: Federico Perini; BSD 3-Clause License,
copyright (c) 2024 Federico Perini. This package uses the
Beliavsky FPM-compatible fork at https://github.com/Beliavsky/fortran-lapack,
pinned to tested commit `9b43beb7644479af46a56e743f61dbc8ae63ac7a`.
The dependency and its original BLAS/LAPACK attribution and license notices
are retrieved separately by FPM, never copied into this package. The manifest
has no system BLAS/LAPACK link flags. The dependency's own portable pure BLAS
and LAPACK chain is used through explicit module interfaces with the same
double-precision real kind.

GPL-2.0-or-later sibling components are used under GPL version 2 for this
combined GPL-2.0-only translation; MIT, Artistic-2.0's GPL-compatible terms,
and the fortran-lapack upstream permissive terms retain their original notices
in their respective dependencies. No license of an upstream component is
replaced by this notice.

## Reference provenance and intentional differences

Actual reference output was generated with R 4.5.0, matrixTests 0.2.3.1, and
matrixStats 1.5.0. No R runtime is used by the Fortran library. Numerical
correctness fixes, resource limits, and omitted R interface behavior are
listed explicitly in README.md and the per-function manifest mappings.
