# NOTICE and provenance

## Upstream R package

This is a modern Fortran translation of computational code from the R package **QZ**, version 0.2-4 (2025-04-12), by Wei-Chen Chen, with LAPACK authors credited by the upstream package.

Upstream metadata is retained in `upstream/DESCRIPTION` and `upstream/NAMESPACE`. The R source files used to establish behavior and coverage are retained in `upstream/R/`.

The upstream package is distributed under the Mozilla Public License 2.0. Its license is retained as `LICENSE`. The LAPACK notice distributed by the upstream package is retained as `LICENSES/LAPACK_LICENSE.txt`.

The translation does not copy the upstream package's C bridge code or fixed-form LAPACK wrapper sources.

## fortran-lapack dependency

Direct LAPACK call chains are supplied through:

- `https://github.com/Beliavsky/fortran-lapack`

That repository is a downstream fork of **Federico Perini's `fortran-lapack`** and uses a modern Fortran implementation of Reference LAPACK. The fork documents its BSD-3-Clause license and its relationship to the original project.

No `fortran-lapack`, BLAS, or LAPACK source is vendored into this package. The dependency is declared in `fpm.toml`.

A commit revision is not pinned in this checkpoint because the sandbox used for translation could inspect the public repository/API but could not perform a Git fetch to resolve and build a specific revision. Pinning a Windows-tested revision is recommended after the first native FPM build in the target repository.

## Translation conventions

The translation keeps the numerical LAPACK operations while replacing R lists/classes, `.Call` bridge code, and S3 dispatch with explicit Fortran derived types and generic interfaces. Printing and dynamic-load hooks are not translated because they are presentation/runtime plumbing rather than computational routines.
