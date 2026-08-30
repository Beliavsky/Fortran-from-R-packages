# NOTICE and provenance

## Upstream package

This work is a modern Fortran translation of computational algorithms from the R
package **vars**, version **1.6-1**.

Upstream metadata identifies:

- Bernhard Pfaff: author and maintainer/creator.
- Matthieu Stigler: contributor.
- Upstream license: `GPL (>= 2)`.
- CRAN publication date in the supplied DESCRIPTION: 2024-03-21.
- The supplied source package declares `NeedsCompilation: no`; its numerical
  implementation is therefore primarily R source rather than an existing
  compiled library.

SHA-256 of the user-supplied `vars-master.zip` used for this translation:

`04390382ca03399ffb44ec7ed72a292ab48c2a9dec18098910587ce71a24467c`

Original package metadata and citation material are retained under `upstream/`:
`DESCRIPTION`, `NAMESPACE`, `CITATION`, and `ChangeLog`.

## Translation provenance

The Fortran implementation follows the numerical definitions and algorithmic
structure in the corresponding upstream R files, while replacing R formula/S3
objects, `lm`/`mlm` dispatch, plotting, printing, and callbacks with typed array
interfaces suitable for Fortran programs.

The translation does not imply endorsement by the upstream authors.

## External shared dependencies

The package uses sibling FPM path dependencies on `rfortran-core` and
`rfortran-linalg`. Those dependencies are not copied, vendored, or embedded here.
No BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, or translated R-package dependency
source is included in this directory.

`rfortran-core` supplies the common `dp` kind and probability-distribution
helpers. `rfortran-linalg` supplies reusable least-squares, inversion,
Cholesky, eigensystem, and SVD operations.

## Upstream citations

The upstream `inst/CITATION` requests citation of:

- Bernhard Pfaff (2008), *VAR, SVAR and SVEC Models: Implementation Within R
  Package vars*, Journal of Statistical Software 27(4).
- Bernhard Pfaff (2008), *Analysis of Integrated and Cointegrated Time Series
  with R*, second edition, Springer, New York.

See the retained `upstream/CITATION` for the complete upstream citation entries.

## License files

`LICENSE` contains GPL-2.0 and `LICENSE-GPL-3.0` contains GPL-3.0. The package
metadata uses the SPDX expression `GPL-2.0-or-later`, corresponding to the
upstream `GPL (>= 2)` declaration.
