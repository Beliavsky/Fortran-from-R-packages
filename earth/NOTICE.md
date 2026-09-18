# NOTICE

This project is a modern Fortran translation of computational portions of the R package **earth**, version 5.3.6.

Upstream package authors listed in `DESCRIPTION`:

- Stephen Milborrow — author and maintainer
- Trevor Hastie — author
- Rob Tibshirani — author

Upstream contributors listed in `DESCRIPTION` include Alan Miller and Thomas Lumley.

The upstream `src/earth.c` records that its MARS code is based on `dmarss` Ratfor from the `mda` package by Trevor Hastie and Rob Tibshirani, with historical R-related modifications attributed there to Kurt Hornik, Friedrich Leisch, and Brian Ripley, followed by earth development by Stephen Milborrow. The upstream source also cites Jerome Friedman's MARS and Fast MARS work.

The upstream package includes `src/leaps.f`, copied from Thomas Lumley's `leaps` package and based on Fortran code by Alan Miller, including Applied Statistics algorithms such as AS 274. This Fortran translation does **not** copy or vendor that dependency-derived source; it uses a package-local modern least-squares implementation. The upstream filename and lineage are documented here for attribution.

The R package declares GPL-3 in `DESCRIPTION`. Individual upstream source files may contain compatible earlier-or-later GPL notices. The maintained translation is distributed under GPL-3; see `COPYING`.

No BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, or translated dependency source is vendored into the maintained Fortran implementation.
