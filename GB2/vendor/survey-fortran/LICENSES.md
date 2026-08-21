# Licensing and provenance

`survey-fortran` is a translation of the computational core of the R package
`survey` 4.5. The upstream package declares `GPL-2 | GPL-3`; translated source
files in `src/` therefore carry `SPDX-License-Identifier: GPL-2.0-only OR
GPL-3.0-only`.

The upstream copyright/provenance notice is preserved verbatim in `LICENSE`.
Full GPL-2.0 and GPL-3.0 license texts are provided in `licenses/`.

## Active vendored dependencies

- `vendor/survival-fortran`: translation of R `survival` 3.8-9, attributed to
  Terry M. Therneau and contributors, licensed `LGPL-2.0-or-later`. Its LGPL
  license texts and a provenance notice are retained inside that directory.
- `vendor/minqa`: supplied Fortran translation of R `minqa`, licensed
  `GPL-2.0-only` according to its bundled metadata.
- `vendor/numDeriv-fortran`: supplied Fortran translation of R `numDeriv`,
  licensed `GPL-2.0-or-later` according to its bundled metadata.

## Reference-only supplied dependencies

- `vendor/MatrixExtra-fortran`: preserved as a supplied reference port. It is
  not in the default FPM dependency graph. Its bundled metadata declares
  `GPL-3.0-only` for the linked MatrixExtra translation.
- `vendor/splines-fortran-reference`: preserved as a supplied standalone
  reference port under `GPL-2.0-or-later`. The active `survival-fortran`
  dependency already contains the compatible translated splines modules, so
  linking this second copy would create duplicate Fortran modules.

No license is relicensed or replaced by this translation. When redistributing
a modified combined work, review the licenses of the components you actually
link and distribute.
