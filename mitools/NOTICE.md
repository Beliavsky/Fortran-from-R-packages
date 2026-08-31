# Notice

This directory is an unofficial modern Fortran translation of computational
code from the R package `mitools` version 2.4.

Upstream metadata identifies:

- Package: `mitools`
- Title: Tools for Multiple Imputation of Missing Data
- Version: 2.4
- Author: Thomas Lumley
- Maintainer: Thomas Lumley
- License field: `GPL-2`
- CRAN publication date: 2019-04-26

The translation preserves the multiple-imputation formulas and numeric
plausible-value/data-container behavior while omitting R-specific S3, formula,
expression-evaluation, database, and interactive infrastructure.

The Fortran source was newly written for this translation on 2026-08-30 and is
distributed under GNU GPL version 2, consistent with the upstream package's
recorded license. No upstream copyright year was asserted in the supplied
metadata, so none is invented here.

`rfortran-core` is a separate sibling FPM dependency under its own MIT license.
Its source is not copied or vendored into this package.

This translation is not endorsed by Thomas Lumley, CRAN, the R Foundation, or
the authors of `rfortran-core`.
