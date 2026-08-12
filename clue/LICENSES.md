# Licenses and provenance

## clue-fortran

This translation is derived from R package `clue` 0.3-68, whose DESCRIPTION
file declares `License: GPL-2`.

The top-level translation is distributed under GPL-2.0-only.  The full GPL v2
text is in `LICENSE`.

Original upstream material is retained under:

```
original/clue-master/
```

## lpSolve-fortran dependency

The vendored dependency under

```
vendor/lpsolve-fortran/
```

is the separately supplied `lpSolve-fortran` 0.1.0 translation and is licensed
LGPL-2.0-only.  Its `LICENSE`, provenance, original sources and notices remain
inside that directory.

The dependency is referenced through FPM as a path dependency and is not
relicensed by this project.
