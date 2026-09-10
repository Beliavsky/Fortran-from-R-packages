# NOTICE and provenance

## Upstream package

This work translates computational code from the R package **signal 1.8-1**
(Date 2024-06-13; upstream package metadata preserved in
`upstream/DESCRIPTION`). The upstream package describes itself as a set of
signal-processing functions translated primarily from Octave and uses the
license declaration `GPL-2`.

Upstream package authors and contributors named in DESCRIPTION include Uwe
Ligges, Tom Short, Paul Kienzle, Sarah Schnackenberg, David Billinghurst,
Hans-Werner Borchers, Andre Carezia, Pascal Dupuis, John W. Eaton, E. Farhi,
Kai Habel, Kurt Hornik, Sebastian Krey, Bill Lash, Friedrich Leisch, Olaf
Mersmann, Paulo Neis, Jaakko Ruohio, Julius O. Smith III, Doug Stewart, and
Andreas Weingessel.

The upstream COPYRIGHTS file records source-level copyrights including Andreas
Weingessel (1995-1997), Friedrich Leisch (1995-1997), Kurt Hornik (1995-1997),
John W. Eaton (1996-1997 and later contributions), Paul Kienzle (1999-2002
periods), Bill Lash (2000), David Billinghurst (2001), Kai Habel (2001-2002),
Andras/Andre Carezia (2002), Julius O. Smith III (2004), and EPRI Solutions,
Inc. (2006), among others. **The authoritative retained notice is
`upstream/COPYRIGHTS`; this summary does not replace it.**

The upstream citation metadata is retained verbatim in `upstream/CITATION`.

## Parks-McClellan / remez provenance

The upstream package includes `src/remez.c`, whose notice attributes the
Parks-McClellan implementation to Jake Janovetz, Copyright (C) 1995, 1998, and
licenses that implementation under the GNU Lesser General Public License,
version 2 or later. The upstream explanatory notice is preserved verbatim as
`upstream/README_Parks-McClellan`.

This Fortran translation **does not copy, vendor, or mechanically translate
`remez.c`**. Its public `remez_filter` routine is an independently written,
deterministic weighted-least-squares approximation with a compatible numerical
purpose, and its coverage entry is explicitly marked partial. The retained
Janovetz notice is included for complete upstream provenance rather than to
claim that the LGPL C implementation is embedded here.

## Other upstream compiled code

The upstream package also contains legacy PCHIP Fortran routines (`dpchim.f`,
`dpchst.f`). They are not copied into this translation. PCHIP behavior is
implemented in new free-form Fortran using a shape-preserving Hermite slope
construction. No legacy fixed-form source is included.

## Dependency provenance

Before translation, the target repository was checked for reusable top-level
packages/shared modules. `rfortran-core`, `rfortran-linalg`, `rfortran-arpack`,
`MASS`, and `pracma` are available there. The translated signal package does not
need to copy or vendor any of them: its required numerical kernels are small and
self-contained, and no external BLAS/LAPACK/ARPACK library is linked.

## Translation license

Maintained Fortran source in this directory carries
`SPDX-License-Identifier: GPL-2.0-only` and is distributed with the GPL version
2 text in `LICENSE`, consistent with the upstream package-level GPL-2
declaration. Retained upstream notice files remain subject to their own stated
terms and attributions.
