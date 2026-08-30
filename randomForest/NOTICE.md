# NOTICE and attribution

This package is a modern Fortran computational translation of the CRAN package `randomForest`, version 4.7-1.2.

Upstream authors listed in `DESCRIPTION`:

- Leo Breiman - author; original Fortran
- Adele Cutler - author; original Fortran
- Andy Liaw - author/maintainer; R port
- Matthew Wiener - author; R port

The upstream native sources contain copyright notices including copyright (C) 2001-2015 Leo Breiman, Adele Cutler and Merck & Co., Inc.; the exact year range varies by source file. `rfsub.f` identifies the original Breiman/Cutler/Merck Fortran and modifications by Andy Liaw and Matt Wiener. `rf.c` identifies Andy Liaw's C driver and Matt Wiener's forest-output modifications.

The original source hashes are recorded in `UPSTREAM.md`. The upstream `DESCRIPTION`, `NAMESPACE`, `CITATION`, `NEWS`, and `COPYING` files are retained under `upstream/`.

The R package declares `License: GPL (>= 2)`. The translated package therefore uses SPDX identifier `GPL-2.0-or-later`. The upstream GPL version 2 text is retained as `LICENSE` and `upstream/COPYING`; GPL version 3 is also included as `LICENSE.GPL-3` because the upstream grant permits later GPL versions.

No claim is made that this translation is authored, endorsed, or supported by the upstream randomForest authors.
