# Licenses and provenance

## relsurv

The source R package `relsurv` 2.3-3 declares `License: GPL (>= 2)`. This
translation is therefore distributed under GPL-2.0-or-later. The complete
upstream package is preserved verbatim under `upstream/relsurv` for attribution
and provenance.

Copyright remains with the upstream authors and contributors, including Maja
Pohar Perme and Damjan Manevski. Translation does not replace or remove upstream
copyright notices.

## survival dependency

The files `src/survival_kinds.f90`, `src/survival_types.f90`,
`src/survival_linalg.f90`, and `src/survival_cox.f90` come from the
user-supplied Fortran translation of `survival` and carry SPDX identifier
`LGPL-2.0-or-later`. The original supplied archive is preserved at
`vendor/survival_f90.zip`.

LGPL-2.0-or-later code is compatible with distribution of this combined project
under the GPL-2.0-or-later terms, while the dependency source retains its own
LGPL notices.
