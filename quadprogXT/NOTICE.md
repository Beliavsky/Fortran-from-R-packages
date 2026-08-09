# Notices

This project is a modern Fortran translation of the computational code in the
R package `quadprogXT` version 0.0.6.

The upstream `quadprogXT` package declares `License: GPL (>= 2)`.  Translated
`quadprogXT`-derived sources are therefore distributed under GPL-2.0-or-later.

The package uses the user-supplied `quadprog-fortran` translation as a vendored
FPM dependency.  That dependency is also GPL-2.0-or-later and retains its own
license/provenance files under `vendor/quadprog-fortran/`.

For strict GNU Fortran builds, the vendored copy of `quadprog_core.f90` has a
warning-only source patch: exact `.EQ.` comparisons against 0 or 1 were
rewritten as exact zero-difference tests such as `abs(x) <= 0`.  This introduces
no numerical tolerance and is intended to preserve the original branch
semantics while avoiding `-Wcompare-reals` under `-Werror`.

The complete supplied upstream R package is preserved under
`original/quadprogXT-master/`.
