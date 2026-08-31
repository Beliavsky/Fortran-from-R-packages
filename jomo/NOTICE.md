# Notices and attribution

## Upstream R package

This work is derived from the computational algorithms in **jomo 2.7-6**,
*Multilevel Joint Modelling Multiple Imputation*.

- Authors: Matteo Quartagno and James Carpenter.
- Upstream package license: GPL-2.
- Upstream package date: 2023-04-13.
- Upstream description cites Carpenter and Kenward (2013),
  DOI 10.1002/9781119942283.
- The package's citation file identifies: Quartagno, M. and Carpenter, J. R.,
  “jomo: A package for Multilevel Joint Modelling Multiple Imputation”.

Copyright in the upstream package remains with its respective copyright
holders.  The modern Fortran translation was prepared in 2026 for the
Fortran-from-R-packages project and is distributed under GPL-2.0-only.

## Native helper provenance retained from upstream

The upstream `src/pdflib.c` and `src/wishart.c` include numerical helper
routines carrying GNU LGPL notices.  Their comments attribute routines to,
among others:

- Barry Brown and James Lovato — original FORTRAN77 random-distribution code.
- Guannan Zhang — original FORTRAN90 chi-square density code.
- John Burkardt — C versions and matrix/Wishart helper routines.

The Fortran files in this translation do **not** copy or vendor those C source
files.  Distribution and dense-linear-algebra kernels were reimplemented in
modern Fortran, while the upstream attribution is retained here and in source
comments because the translated algorithms follow the same mathematical
operations.  `LICENSE-LGPL-2.1` is included to preserve the LGPL notice context
conservatively; the package as a whole is distributed under GPL-2.0-only.

## No vendored dependencies

No source from `lme4`, `survival`, `MASS`, `ordinal`, `tibble`, BLAS, LAPACK,
ARPACK, `rfortran-core`, `rfortran-linalg`, or other translated packages is
included in this directory.
