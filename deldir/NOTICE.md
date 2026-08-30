# Notices and provenance

This project is a modern Fortran translation of the computational code in the R package `deldir`, version 2.0-4 (2024-02-27), authored and maintained upstream by Rolf Turner.

The upstream DESCRIPTION declares:

- Package: deldir
- Version: 2.0-4
- License: GPL (>= 2)

The Delaunay/Dirichlet implementation originates in Rolf Turner's implementation of the second iterative algorithm of Lee and Schacter, *International Journal of Computer and Information Sciences*, Vol. 9, No. 3 (1980), pp. 219-242. Upstream source comments and documentation contain additional historical attribution.

The original package source used for this translation is retained under `upstream/deldir-2.0-4/`. The maintained Fortran kernel in `src/deldir_kernel.f90` is a modernization of those upstream Fortran routines and therefore remains covered by the upstream GPL terms.

Changes in this translation include free-form module organization, explicit typing and interfaces, a package-wide `real64` kind, native Fortran error handling, typed result objects, native wrappers for R-side computational orchestration, half-plane construction of rectangular-window tile polygons, and FPM build/test support.

Version 0.2.0 adds the computational behavior of upstream `tile.list(..., clipp=...)`, `doClip()`, and `tileInfo(..., clipp=...)` through the separately translated `polyclip-fortran` v0.1.0 package. That dependency is vendored under `vendor/polyclip-fortran/`, including its BSL-1.0 `LICENSE`, notices, translation notes, and retained upstream `polyclip` / Clipper 6.4.0 source. The deldir project remains distributed under GPL-2.0-or-later; the dependency retains its own license.

See `COPYING` and `vendor/polyclip-fortran/LICENSE`.
