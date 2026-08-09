# Upstream provenance

Translated from the supplied `limSolve-master` package, version 2.0.3.

Upstream authors include Karline Soetaert, Karel Van den Meersche, and Dick van
Oevelen.  The package also credits Charles L. Lawson and Richard J. Hanson for
`inverse.f`, Jack Dongarra for `solve.f`/`inverse.f`, and Cleve Moler for
`solve.f`.

The upstream DESCRIPTION declares `License: GPL` without a version qualifier;
this release preserves that declaration rather than inventing a GPL version.
The upstream `inst/COPYRIGHTS` notice for bundled LAPACK-derived material is
preserved in the complete original source tree.

Vendored dependencies preserve their own licenses:

- `lpSolve-fortran`: LGPL-2.0-only
- `quadprog-fortran`: GPL-2.0-or-later

The vendored quadprog source contains a warning-only portability patch replacing
exact `.EQ.` real comparisons with exact zero-difference inequality tests so it
can compile under `-Werror=compare-reals`.  No numerical tolerance was added.
