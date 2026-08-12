# Licenses and provenance

## adagio-fortran

The translated `adagio` code in `src/`, tests/examples written for this port,
and associated documentation are distributed under **GPL-3.0-or-later**, in
accordance with the upstream R package declaration `License: GPL (>= 3)`.

The full GNU GPL version 3 text is in `LICENSE`.

Upstream attribution:

- Package: adagio 0.9.2
- Author/Maintainer: Hans W. Borchers
- Original package source retained at `original/adagio-master/`

## Vendored lpSolve-fortran

`vendor/lpSolve-fortran-v0.1.0/` is the user-supplied modern Fortran
translation of the R package `lpSolve` / lp_solve computational interface.  It
is a separate FPM dependency and declares **LGPL-2.0-only**.

Its own `LICENSE`, `NOTICE.md`, `UPSTREAM_PROVENANCE.md`, and source-level SPDX
headers are retained unchanged.  The adagio translation does not relicense
those dependency sources.

The final source distribution therefore contains GPL-3.0-or-later application
code linked against an LGPL-2.0-only library dependency, with the two bodies of
source and their notices kept separately identifiable.
