# Upstream provenance and licenses

## survival

Source archive: `survival-master.zip`

- Package: survival
- Version: 3.8-9
- Date: 2026-07-07
- Authors: Terry M. Therneau; Thomas Lumley; Elizabeth Atkinson; Cynthia Crowson
- Upstream license declaration: `LGPL (>= 2)`

The complete supplied source tree is retained under
`original/survival-master/`. The archive did not include a standalone COPYING
or LICENSE text, so this translation preserves the package DESCRIPTION license
declaration rather than inventing a missing upstream file.

## splines-fortran

The user supplied `splines-fortran-v0.1.0`, a modern Fortran translation of R
`splines 2.0-7`. It is retained under `vendor/splines-fortran-v0.1.0/` and under
`original/splines-fortran-v0.1.0/` for provenance. Its translated sources are
GPL-2.0-or-later and its full GPL-2 license text is preserved inside that
package.

The survival-derived `.f90` files remain LGPL-2.0-or-later. The root FPM target
links to the GPL spline dependency for `pspline_basis`, so redistribution of
the linked combined work is governed by GPL-compatible terms; the root
`fpm.toml` therefore declares GPL-2.0-or-later.
