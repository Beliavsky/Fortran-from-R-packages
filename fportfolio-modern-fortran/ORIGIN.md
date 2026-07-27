# Origin and licensing

This project is a clean modern Fortran translation of computational ideas and
formulas in the R package `fPortfolio`, version 4023.84, packaged 2023-04-23.

Original metadata is retained in `original_metadata/`:

- `DESCRIPTION`
- `NAMESPACE`
- `ChangeLog`

The original package declares `License: GPL (>= 2)`. The translation therefore
uses GPL-2.0-or-later. Every Fortran source file includes the corresponding SPDX
identifier and an explicit GPL version 2-or-later notice.

The original package delegates several algorithms to separately licensed
solver and statistics packages. Their source code is not copied here. The
Fortran project supplies self-contained numerical algorithms for the translated
portfolio workflows and documents where those algorithms are analogues rather
than reproductions of external backends.
