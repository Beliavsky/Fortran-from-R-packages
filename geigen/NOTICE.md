# Notices and provenance

This directory is a modern Fortran translation of the computational code in the R package `geigen` version 2.4 (2026-07-20), authored by Berend Hasselman and maintained by Jonathan M. Lees, with LAPACK authors credited by the upstream package.

The upstream package is licensed under GPL version 2 or later. The GPL-2 text is retained in `LICENSE`. The upstream copyright notice for LAPACK-derived source files is retained verbatim in `NOTICE_UPSTREAM.txt`, and the LAPACK license is retained in `LICENSES/LAPACK_LICENSE.txt`.

The upstream R package distributed fixed-form Fortran wrappers and selected LAPACK sources. Those sources are **not** copied into this translation. Generalized eigenvalue, generalized Schur, SVD, eigenvalue, and RQ operations use the FPM-compatible Beliavsky fork of `fortran-lapack`:

`https://github.com/Beliavsky/fortran-lapack`

`fortran-lapack` was originally authored by Federico Perini. The Beliavsky fork is used because it is arranged for the target FPM workflow. LAPACK remains the work of its respective authors and copyright holders.

The public `la_lapack` module in the dependency does not expose a `ggsvd3` generic. Consequently, this translation does not vendor `D/ZGGSVD3` or `D/ZGGSVP3`. The GSVD is reconstructed from public LAPACK primitives exposed by the dependency (`GESVD`, `SYEV`/`HEEV`, `GERQF`, and `ORMRQ`/`UNMRQ`) and is tested through its defining decomposition identities.

Upstream R source files used to establish behavior and coverage are retained under `upstream/R/` solely for provenance and traceability.
