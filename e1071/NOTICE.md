# Notice and provenance

This project is a modern Fortran translation of computational functionality from the R package `e1071` version 1.7-17.

The upstream package authors listed in `DESCRIPTION` are David Meyer, Evgenia Dimitriadou, Kurt Hornik, Andreas Weingessel and Friedrich Leisch. The upstream package also credits Chih-Chung Chang and Chih-Chen Lin for the bundled LIBSVM C++ code.

The upstream package declares `GPL-2 | GPL-3`. Accordingly, this translation is distributed under GPL-2.0-only OR GPL-3.0-only. Copies of GPL version 2 and GPL version 3 are included under `licenses/`.

The unmodified upstream package snapshot used for the translation is retained under `upstream/` for provenance, license review and parity auditing.

The translated `proxy-fortran` dependency is vendored under `dependencies/proxy-fortran`; its own notice, licenses and retained upstream sources remain with that dependency.

The modern Fortran implementation is a translation/reimplementation and is not presented as an official release of the upstream `e1071` authors or the LIBSVM authors.
