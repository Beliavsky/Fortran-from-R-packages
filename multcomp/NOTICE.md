# Notices and provenance

This project is a modern Fortran translation of the computational portions of
`multcomp` 1.4-32, by Torsten Hothorn, Frank Bretz, Peter Westfall, and
contributors. The upstream R package is distributed under GPL-2. The complete
source snapshot used for this translation is retained under `upstream/`.

The translated simultaneous-inference routines use the previously translated
`mvtnorm-fortran` 1.4.2 package under `dependencies/mvtnorm-fortran/`.
`mvtnorm` is also distributed under GPL-2. Its original source and provenance
files are retained inside that dependency.

The Fortran translation is a derivative work and is distributed under
GPL-2.0-only. See `LICENSE`.
