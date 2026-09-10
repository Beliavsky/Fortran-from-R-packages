# Notices and provenance

Upstream package: `svd` 0.5.8, "Interfaces to Various State-of-Art SVD and Eigensolvers".
Upstream project: https://github.com/asl/svd

Upstream credited authors/copyright holders include Anton Korobeynikov, Rasmus Munk
Larsen / Stanford University, and Lawrence Berkeley National Laboratory. Exact upstream
metadata is retained under `upstream/`, and the package-wide extracted-source checksums
are in `provenance/UPSTREAM_SHA256SUMS.txt`.

The upstream DESCRIPTION declares `BSD_3_clause + file LICENSE`. The computational R/C
wrapper files also contain explicit GNU GPL version 2-or-later notices. To avoid weakening
those file-level notices, maintained Fortran source in this translation is distributed
under GPL-2.0-or-later; the GPL text is in `LICENSE`. The exact upstream CRAN-style
`LICENSE` metadata remains at `upstream/LICENSE`.

The upstream source bundle contains PROPACK (copyright 2005 Rasmus Munk Larsen, Stanford
University) and nuTRLan (copyright 2009 Lawrence Berkeley National Laboratory). Neither
library is copied into, compiled by, or linked into this translation. Their exact license
texts and the upstream COPYRIGHTS file are retained under `provenance/` for attribution
and historical provenance only.

Numerical dependencies are sibling FPM path dependencies:

- `../RSpectra` for real dense/matrix-free SVD and symmetric eigenvalue computation.
- `../rfortran-linalg` for conventional complex dense SVD.

Those dependencies retain their own licenses and are not vendored here.
