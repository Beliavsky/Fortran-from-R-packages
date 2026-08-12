# Licensing and provenance

## Upstream ppso package

The supplied upstream package is ppso 0.9-99994 by Till Francke.
Its `DESCRIPTION` file states exactly:

```text
License:  Unlimited
```

The supplied repository does not contain a separate LICENSE or COPYING file.
The complete upstream package is retained under `original/`, including the
original `DESCRIPTION`, R sources, documentation, and README.

## Fortran translation

This translation preserves the upstream `Unlimited` license declaration and
does not impose a more restrictive license on translated computational code.
Where comments identify an algorithm or source routine, those comments are for
provenance and traceability.

No Rmpi, lhs, or rgl source code has been copied into the Fortran library.
Latin-hypercube initialization and the standalone RNG are new self-contained
Fortran implementations.
