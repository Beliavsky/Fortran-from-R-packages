# Licensing

## sadists-derived Fortran code

The supplied R package `sadists` 0.2.6 identifies its license as LGPL-3, and
its R source headers grant redistribution/modification under GNU LGPL version
3 or, at the recipient's option, any later version.

The Fortran translation derived from those files is therefore distributed as:

**SPDX-License-Identifier: LGPL-3.0-or-later**

The full GNU LGPL version 3 text is in `LICENSES/LGPL-3.0.txt`.

## PDQutils-fortran dependency

Version 0.2.0 depends on the standalone `pdqutils-fortran` 0.1.0 translation,
vendored under `vendor/pdqutils-fortran/`. It is also distributed under
LGPL-3.0-or-later. Its own license, notices, API map, and porting notes are
retained inside that directory.

The compiled sadists approximation layer no longer contains a private copy of
PDQutils Edgeworth, Cornish-Fisher, AS269, normal-probability, or
moment/cumulant algorithms.

## Other mathematical algorithms

Sadists-specific gamma, Poisson, noncentral-chi-square, special-function, and
Poisson-mixture implementations were written for this port from standard
mathematical identities/algorithms; they do not vendor R, hypergeo, or
orthopolynom code.
