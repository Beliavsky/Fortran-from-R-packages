# Licensing and provenance

## mstate

The upstream `mstate` package declares `GPL (>= 2)` and is distributed under
the GNU General Public License version 2 or, at the user's option, any later
version. This Fortran translation is distributed under the same
GPL-2.0-or-later terms.

The complete upstream package is retained in `upstream/mstate-0.3.3/`.
Algorithm provenance is especially direct for `src/mstate_msfit.f90`, which is
a Fortran translation of upstream `src/agmssurv.c`.

## survival dependency

The user-supplied `survival_f90.zip` source bundle is retained verbatim under
`vendor/survival_f90/`. Its Fortran sources declare `LGPL-2.0-or-later` in
their SPDX headers. The Cox-related modules copied into `src/` retain those
headers unchanged. LGPL-2.0-or-later code can be combined with this
GPL-2.0-or-later work while retaining the LGPL notices on the survival files.

## relsurv dependency

`relsurv` 2.3-3 declares `GPL (>= 2)`. The complete translated dependency is
retained as `vendor/relsurv-fortran-v0.2.0.zip`. The modules copied into `src/`
(`relsurv_kinds`, `relsurv_ratetable`, and `relsurv_parsers`) remain under
GPL-2.0-or-later and are license-compatible with mstate.

No upstream copyright notices or license declarations have been removed from
the retained source/dependency archives.
