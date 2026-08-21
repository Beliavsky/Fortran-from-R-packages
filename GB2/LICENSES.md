# Licensing and provenance

## GB2 translation

`GB2-fortran` translates the computational core of R package `GB2` 2.1.2
(2025-08-17), authored by Monique Graf and Desislava Nedyalkova. Upstream
`DESCRIPTION` declares `License: GPL (>= 2)`. Translated files under `src/`
carry `SPDX-License-Identifier: GPL-2.0-or-later`.

The unmodified upstream package snapshot used for translation is retained at
`upstream/GB2-2.1.2/`. Full GPL-2.0 and GPL-3.0 texts are in `licenses/`.

## Active dependency

The default FPM graph links `vendor/survey-fortran`, the previously supplied
translation of R `survey`. Its own `LICENSES.md` records the upstream and
vendored dependency licenses. The survey translation declares
`GPL-2.0-only OR GPL-3.0-only`; for a GPL-2.0-or-later GB2 combined work, its
GPL-2.0 option is compatible.

## Supplied reference-only dependency translations

The user supplied three additional translations. They are retained for
provenance/reference, but are deliberately **not** in this package's FPM
link dependency graph:

- `vendor/cubature-fortran-reference`: `GPL-3.0-or-later` according to its
  bundled `fpm.toml`.
- `vendor/hypergeo-fortran-reference`: `GPL-2.0-only` according to its bundled
  `fpm.toml`.
- `vendor/numDeriv-fortran-reference`: `GPL-2.0-or-later` according to its
  bundled `fpm.toml`.

The supplied cubature and hypergeo translations cannot both be linked into a
single binary under their stated licenses: GPL-3.0-or-later and GPL-2.0-only
have no common license version. The GB2 translation therefore implements the
small numerical pieces it needs internally (one-dimensional adaptive
Gauss-Kronrod integration and the convergent real 3F2(1) series). This avoids
creating a conflicted binary license graph while preserving the supplied
ports separately as an aggregation/reference snapshot.

The active survey dependency already contains the numerical-derivative support
needed by its own routines. GB2's direct derivatives are analytic or use local
finite-difference/Jacobian helpers, so the separately supplied numDeriv port is
not needed in the default link graph.
