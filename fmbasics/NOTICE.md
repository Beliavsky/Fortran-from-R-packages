# Notice

This project is a modern Fortran translation of the computational portions of
`fmbasics` 0.3.99 by Imanuel Costigan and contributors.

The upstream package declares the GPL-2 license in its DESCRIPTION metadata.
This translation is distributed under GPL-2.0-only. The complete license text
is in `LICENSE`.

Original R source files, tests, DESCRIPTION, NAMESPACE, NEWS and README are
retained under `original/`. Their presence is for attribution, provenance and
comparison; the translated Fortran sources are the files under `src/`, `test/`,
`example/` and `app/`.

The internal calendar and CDS bootstrap implementations are adaptations for a
self-contained Fortran package and should not be presented as the original
`fmdates` or `credule` implementations.
