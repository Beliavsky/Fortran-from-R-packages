# Licenses

## linprog translation

The attached R package `linprog` 0.9-6 declares `License: GPL (>= 2)`.
The modern Fortran translation in `src/`, `test/`, and `example/` is therefore
provided under **GPL-2.0-or-later**.  The original package sources are retained
under `original/linprog-master/` for provenance.

## lpSolve-fortran dependency

The user-supplied `lpSolve-fortran-v0.1.0` translation is vendored unchanged
under `vendor/lpsolve-fortran/`.  It declares **LGPL-2.0-only** and retains its
own `LICENSE`, `NOTICE.md`, and provenance documentation.  It is an FPM path
dependency and is not copied into the linprog module namespace.
