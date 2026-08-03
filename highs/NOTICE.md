# Notice

This project is a modern Fortran translation of the computational interface exposed by the R package `highs` 1.14.0-2.

The R package is licensed under GPL version 2 or later. The Fortran frontend, dynamic loader, bridge, tests, examples, and documentation in this project are distributed under GPL-2.0-or-later.

The numerical optimizer is HiGHS 1.14.0, retained under `upstream/highs-r/inst/HiGHS` and built as an external runtime backend. HiGHS is distributed under the MIT License; its original `LICENSE.txt` is retained unchanged. The bridge links against HiGHS but does not change the solver algorithms.

The complete attached R-package snapshot is retained in `upstream/highs-r` for provenance. Files in that directory remain under their original licenses.
