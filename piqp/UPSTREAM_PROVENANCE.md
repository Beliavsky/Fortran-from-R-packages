# Upstream provenance

* Upstream R package: `piqp` 0.6.2.
* Upstream package description: R interface to the PIQP proximal interior-point
  quadratic-programming solver.
* PIQP-R license: BSD 2-Clause (`BSD_2_clause + file LICENSE`).
* The bundled C++ PIQP headers contain BSD 2-Clause notices and copyrights for
  EPFL/INRIA contributors.
* Original uploaded package is preserved verbatim in `original/piqp-master/`.
* User-supplied Matrix Fortran translation is preserved in
  `vendor/Matrix-fortran/` and retains its GPL-3.0-only license and provenance.

The translated solver modules in `src/` retain BSD-2-Clause SPDX identifiers.
The linked distribution is GPL-3.0-only because of the Matrix dependency; see
`NOTICE.md`.
