# Licensing and provenance

## genalg-fortran

This project is a translation and derivative work of the computational code
in the R package `genalg` 0.2.1.

The original package DESCRIPTION declares:

```text
Package: genalg
Version: 0.2.1
Author: Egon Willighagen and Michel Ballings
License: GPL-2
```

Accordingly, the translated source is distributed under the GNU General
Public License version 2. The full GPL v2 text is in `LICENSE`.

The original package material used for the translation is retained under
`original/`, including `DESCRIPTION`, `NAMESPACE`, `ChangeLog`, the R sources,
and manual pages.

No third-party numerical library code is incorporated into the translated
Fortran source. The standalone RNG in `src/genalg_rng.f90` is a direct small
implementation of the classic Park-Miller minimal-standard recurrence and is
written as part of this translation.
