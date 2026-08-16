# Licenses

## mbbefd-fortran translated source

The supplied R package `mbbefd` 0.8.14 declares `License: GPL-2`. The translated
mbbefd-derived source in `src/` is distributed under **GPL-2.0-only** and carries
SPDX identifiers accordingly. The complete upstream package is retained under
`upstream/mbbefd-master/`.

A copy of GNU GPL version 2 is in `LICENSE`.

## Vendored dependencies

The user supplied the following pre-existing Fortran translations, retained
under `vendor/` with their own notices and license files:

- `fitdistrplus-fortran` 0.1.0 — GPL-2.0-or-later
- `actuar` 0.1.0 — GPL-2.0-or-later
- `alabama` 0.1.1 — GPL-2.0-or-later
- nested dependencies of `alabama`, including the supplied `numDeriv` and
  `roptim` ports, under the licenses recorded in that vendor tree

Those files remain governed by their respective notices. The combined package
is distributed under terms compatible with GPL version 2.
