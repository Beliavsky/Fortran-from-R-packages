# Licensing and provenance

This project is a modern Fortran translation of the computational code in the
R package **msm 1.8.2**, by Christopher Jackson.  The source package declares
`License: GPL (>= 2)`, so this translation is distributed as
**GPL-2.0-or-later**.

The original `DESCRIPTION`, `NAMESPACE`, native C headers/sources, and selected
computational R sources are retained under `original/` for provenance and
algorithm auditing.  They are not compiled as part of the FPM library.

The matrix-exponential implementation in `msm_linalg.f90` is a fresh Fortran
implementation of the scaling-and-squaring / Pade approach used for the same
purpose by `msm` and its `expm` dependency; no third-party source is linked or
vendored beyond BLAS/LAPACK.
