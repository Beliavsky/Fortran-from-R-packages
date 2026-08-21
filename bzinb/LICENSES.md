# Licensing and provenance

`upstream/` is a verbatim copy of R package `bzinb` 1.0.8 as supplied by the
user. Its DESCRIPTION declares `License: GPL-2`.

The Fortran translation follows those probability laws and algorithms and is
released under GPL-2.0-only. The GPL version 2 text is included as `LICENSE`.

In v0.2.0, `src/bzinb_em.f90` is a direct numerical translation of the
GPL-licensed upstream `src/expt.cpp`, `src/opt.cpp`, and `src/em.cpp`, with Rcpp
and Boost replaced by native Fortran arrays and the package's native
digamma/trigamma implementations.

The translation does not link Rcpp or Boost and has no mandatory external
runtime library dependency.
