# rfortran-arpack

This internal FPM package supplies the double-precision symmetric and
nonsymmetric ARPACK routines used by translated packages. It avoids a system
ARPACK, BLAS, or LAPACK installation by using the repository's pinned
pure-Fortran LAPACK backend.

The double-precision ARPACK-NG 3.9.1 sources were taken from the official
release and mechanically converted from fixed source form to free source form.
The `rfortran_arpack` module supplies explicit interfaces for `dsaupd`,
`dseupd`, `dnaupd`, and `dneupd`.

See `LICENSE` for the upstream BSD license and attribution.
