# NOTICE

This directory is a modern Fortran translation of computational code from the
R package `kde1d` version 1.2.2 by Thomas Nagler and Thibault Vatter.

The supplied upstream package is licensed under the MIT License. Its original
R license metadata is preserved in `LICENSE` and `upstream/LICENSE`. The full
MIT license text is provided in `LICENSE-MIT`. The bundled upstream
`kde1d-cpp` license is preserved verbatim at
`upstream/include/kde1d-cpp/LICENSE`.

Translated numerical code derives from the R and C++ sources in the supplied
archive. See `PROVENANCE.md` for numerical provenance and documented
implementation differences.

The package uses the sibling `rfortran-core` package from
`Fortran-from-R-packages` as an external FPM dependency. No source from that
dependency is copied into this package.
