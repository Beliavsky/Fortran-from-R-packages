# Dependency notes

## Linked source

### tsdistributions Fortran translation

The distribution, random-generation, moment, profiling, optimization, and
special-function modules from `tsdistributions-fortran-v0.1.0` are compiled
into this package. Those translated sources carry GPL-2.0-only notices.

### ghyp Fortran translation

The generalized-hyperbolic numerical modules carried by the distribution layer
retain GPL-2.0-or-later file headers. GPL-2.0-or-later code may be distributed in
this GPL-2.0-only combined work under GPL version 2.

## Provenance-only source

### nloptr Fortran translation

`provenance/dependencies/nloptr-fortran(2).zip` is the exact attachment supplied
for comparison. It is not listed in `fpm.toml`, the Makefile, any shell script,
or any Fortran `use` statement.

The archive declares LGPL-3.0-or-later. LGPL version 3 incorporates GPL version
3 terms, which cannot be satisfied by a combined work whose upstream license is
GPL-2.0-only. Keeping the archive unlinked as separate provenance avoids that
license conflict.

The interface review found that the supplied optimizer could conceptually cover
bounded and constrained minimization, but replacing the current optimizer would
not change the statistical API. A future upstream relicensing to GPL-2.0-or-
later or GPL-3.0-or-later would permit reevaluating that integration.
