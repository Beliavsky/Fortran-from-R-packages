# Upstream provenance

This project translates the computational portions of the R package `neldermead`
version 1.0-13 (2026-01-25), supplied in `neldermead-master.zip`.

The R package is itself a port of Michael Baudin's Scilab Nelder-Mead component.
The upstream DESCRIPTION declares `License: CeCILL-2` and dependencies on
`optimbase (>= 1.0-9)` and `optimsimplex (>= 1.0-7)`.

The supplied upstream archive is copied verbatim under `original/neldermead-master/`.
It references a `COPYING` file in source headers, but that file was not present in
the supplied archive. `LICENSE` therefore preserves the declared license and
points to the authoritative CeCILL-2 English license rather than inventing text.
