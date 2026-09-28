# NOTICE and provenance

This directory is a modern free-form Fortran translation of computational behavior exposed by the R package **rmumps** version 5.2.1-43.

Upstream rmumps wrapper:
- Author/maintainer: Serguei Sokol, with contributors listed in `upstream/DESCRIPTION`.
- Wrapper copyright: 2015-2026 INRAE/INSA/CNRS.
- Upstream wrapper license: GPL (>= 2).
- Upstream project: https://github.com/sgsokol/rmumps/
- MUMPS project cited by upstream: http://www.mumps-solver.org/

The original R package bundles MUMPS and ordering-library sources. Those third-party sources are **not** copied or vendored in this translation. Their upstream copyright and license notices are preserved in `UPSTREAM_COPYRIGHTS` for provenance.

Current translation status: the maintained Fortran API uses a package-local dense LU factorization with partial pivoting. It preserves matrix construction from triplets, cached numerical factorization, dense vector and matrix solves, sparse-triplet RHS solves, transposed solves using the same factors, inverse and inverse-transpose operations, determinant operations, triplet access, value updates, symmetry modes, dimensions, permutation-code state behavior, and the public ICNTL/CNTL/KEEP state-management methods. MUMPS-shaped information arrays are also exposed with the status and determinant information available from the fallback backend. Public wrapper state for `copy`, `ncore`, `rhs`, and `mrhs` is represented, and user permutation requests can be validated and stored with the same ICNTL state changes and factor invalidation performed by upstream rmumps.

The implementation does **not** reproduce MUMPS sparse multifrontal factorization, ordering algorithms, application of user-supplied MUMPS permutations, numerical effects of most MUMPS controls, MUMPS-only diagnostic fields, parallel or distributed execution, sparse memory scaling, or Rcpp aliasing/object semantics. Its dense backend is intended as a portable computational translation/fallback and validation target, not as a performance substitute for MUMPS on large sparse systems.

No MUMPS, BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, or translated dependency source is vendored in this package. The current slice also does not require a system BLAS/LAPACK installation.

No source from Federico Perini's `fortran-lapack` or the Beliavsky FPM-compatible fork is copied into this package, and this current slice does not depend on it. If a later revision replaces or supplements the package-local fallback with direct LAPACK calls, the dependency and attribution required by the project instructions must be added explicitly.
