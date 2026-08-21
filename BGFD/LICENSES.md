# Licensing and attribution

## BGFD

The upstream R package BGFD 0.1 declares `License: GPL (>= 2)`. The translated
BGFD code in this repository is therefore distributed under
**GPL-2.0-or-later**. See `LICENSE`.

Original authors:

- Michail Tsagris
- Muhammad Imran
- M.H. Tahir

The supplied upstream package snapshot is retained in `upstream/BGFD-0.1/`.

## AdequacyModel Fortran translation

BGFD uses `AdequacyModel::goodness.fit` for its `m*` fitting routines. The
previously translated AdequacyModel numerical modules are vendored directly in
`src/` so this FPM project has no external package-manager dependency. Those
files retain their GPL-2.0-or-later SPDX headers and attribution.

The accompanying notices from that port are retained in
`third_party/AdequacyModel-fortran-v0.1.0/`.

Because both BGFD and the vendored AdequacyModel translation permit GPL version
2 or later, the combined source tree is license-compatible as
GPL-2.0-or-later.
