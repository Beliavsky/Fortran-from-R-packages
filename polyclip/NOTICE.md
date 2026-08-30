# Notices and provenance

This project is a modern Fortran translation of the computational geometry
provided by the R package **polyclip 1.10-7**.

The R package is an R port of Angus Johnson's **Clipper** library and records
that it was built from Clipper C++ version 6.4.0. The retained upstream source
in `upstream/` is included for licensing, provenance, and parity review.

Original Clipper code:

- Author: Angus Johnson
- Copyright: Angus Johnson 2010-2015
- License: Boost Software License 1.0

R polyclip port/contributors include Adrian Baddeley, Kurt Hornik, Brian D.
Ripley, Elliott Sales de Andrade, Paul Murrell, Ege Rubak, and Mark Padgham.
See `upstream/DESCRIPTION` for the complete package metadata.

The maintained Fortran source is a translation/derived work distributed under
the same Boost Software License 1.0. See `LICENSE`.
