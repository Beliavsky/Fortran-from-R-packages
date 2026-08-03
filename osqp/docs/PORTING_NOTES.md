# Porting notes

## Why a backend bridge is used

The attached R package is a language binding to the OSQP C solver. Replacing OSQP's ADMM, scaling, QDLDL factorization, polishing, certificates, adaptive rho logic, and stopping rules with a small independent Fortran optimizer would not be an equivalent translation. This project therefore provides an idiomatic Fortran API over the same bundled OSQP 1.0.0 engine.

## Runtime loading

The frontend calls a small C dynamic loader compiled by FPM. It searches:

- `OSQP_FORTRAN_BRIDGE`
- the current directory
- `backend/bin` on Windows
- `backend/lib` on Unix

Consequently, `fpm build` does not require a prebuilt native library or linker search path. A missing backend is represented by `osqp_backend_unavailable` and is reported without a link failure.

## Sparse representation

Public CSC arrays are one-based and sorted within columns. The backend bridge receives temporary zero-based arrays. The bridge copies all setup data before calling OSQP, so Fortran temporaries do not outlive their storage.

Only the upper triangle of `P` is retained. Matrix value updates preserve the original sparsity pattern, as required by OSQP.

## Bounds

IEEE infinities and values beyond OSQP's finite sentinel are clamped to `+/-OSQP_INFTY` inside the bridge. The Fortran model retains the original values.

## Settings updates

OSQP permits some settings only during setup. `osqp_update_settings` delegates to the native OSQP update functions; setup-only fields are therefore ignored by the backend exactly as documented by OSQP. The Fortran solver then reads back the active settings.

## Backend configuration

The offline backend build uses:

- double precision
- 32-bit integer indices, matching the R package
- built-in QDLDL direct linear solver
- printing, timing, and interrupt support
- matrix updates enabled
- code generation and derivative support disabled

The generic printing header replaces the R-specific `Rprintf` header copied by the R package's configure script.
