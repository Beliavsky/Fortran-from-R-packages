# Licensing notice

This distribution combines two separately licensed code bases.

* The translated PIQP / PIQP-R derived source in `src/` is distributed under
  the BSD 2-Clause License. The upstream R package declares
  `BSD_2_clause + file LICENSE`, and the bundled PIQP C++ headers carry BSD
  2-Clause notices.
* `vendor/Matrix-fortran/` is the user-supplied Fortran translation of the R
  Matrix package and is GPL-3.0-only.

Because the default FPM target links `matrix-fortran` to provide the sparse CSC
adapter, redistribution of the combined linked target is under GPL-3.0-only
terms. The BSD license remains applicable to the PIQP-derived source files
individually.

The full upstream PIQP-R package is retained in `original/piqp-master/`.
