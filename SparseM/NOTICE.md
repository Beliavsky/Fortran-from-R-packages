# NOTICE and provenance

This work translates computational portions of **SparseM 1.84-2** from R and legacy Fortran to modern free-form Fortran.

Upstream package metadata identifies:

- Roger Koenker - creator/author
- Pin Tian Ng - contributions to sparse QR code
- Yousef Saad - author of SPARSKIT2 code used upstream
- Ben Shaby - author of upstream `chol2csr`
- Martin Maechler - Cholesky tweaks and S4 work

Upstream license: **GPL (>= 2)**. This translation is distributed as **GPL-2.0-or-later**. The full GPL version 2 text is in `LICENSE`.

Upstream source package snapshot supplied by the user: `SparseM-master.zip`, package version 1.84-2, with DESCRIPTION publication date 2024-07-17. The files `upstream/DESCRIPTION`, `upstream/NAMESPACE`, and `upstream/R/SparseM.R` are retained for attribution and coverage traceability. The upstream `inst/extdata/rua_32_ax.rua` Harwell-Boeing example is retained under `upstream/inst/extdata/` and copied to `test/data/` solely as a deterministic parser fixture.

The original SparseM package contains fixed-form Fortran derived in part from the Ng-Peyton sparse Cholesky implementation and SPARSKIT2. The maintained source does **not** vendor those fixed-form files. `src/sparsem_ng_order.f90` is a modern free-form translation of the multiple-minimum-degree ordering stage (`ORDMMD`, `GENMMD`, `MMDELM`, `MMDINT`, `MMDNUM`, and `MMDUPD`). `src/sparsem_ng_symbolic.f90` translates the central SFINIT/SYMFCT symbolic path, including elimination-tree construction/postordering, column counting, child reordering, fundamental-supernode discovery, and compressed supernodal row-index construction. `src/sparsem_ng_numeric.f90` translates the INPNV/BFINIT preparation, Ng-Peyton left-looking block numerical factorization, cache splitting, the level-1 update path plus paired-target level-2/4/8 update scheduling, and the BLKSLF/BLKSLB/BLKSLV triangular-solve algorithms. `chol_statistics` reproduces the computational statistics calculated by the upstream `LSTATS` helper without its printing side effects.

These translated routine families preserve upstream algorithm/source credits to Joseph W. H. Liu, Esmond G. Ng, Barry W. Peyton, Guodong Zhang, the University of Waterloo SPARSPAK-A release (January 1984), and the Oak Ridge National Laboratory Mathematical Sciences Section. The original fixed-form sources were used only as validation oracles during translation; they are not included in this package. Deterministic comparison tests matched the modern ordering and symbolic outputs against the original routines on 40 graph patterns. End-to-end development-oracle tests using the original SparseM `chol` driver matched permutation vectors, symbolic structures, cache splits, and final `LNZ` values exactly at levels 1, 4, and 8 for tridiagonal, star, grid, and irregular sparse positive-definite cases. Direct `FNSPLT` and `FNTSIZ` oracle tests also matched cache partitioning and temporary-workspace sizing, including a 200-column supernode split across multiple cache blocks. All four modern numerical levels were additionally checked for equivalent factors/solutions over a larger deterministic matrix set, including a package test spanning twelve additional sparse positive-definite patterns. The supplied legacy level-2 path produced divergent factors on some of those patterns, so it is not treated as a bit-for-bit compatibility oracle; the maintained level-2 path implements the same intended numerical operation as levels 1, 4, and 8. A compact upstream-reference grid case is retained as `test/test_ng_reference.f90` so routine FPM testing can enforce ordering, symbolic-layout, cache/workspace, and numerical-factor invariants without copying the legacy fixed-form sources.

The public vector and matrix right-hand-side APIs now call the translated supernodal triangular-solve routines directly when an Ng-Peyton factor is available; the CSR mirror remains for interoperability and inspection rather than as the primary solve path.

No source from `rfortran-core`, `rfortran-linalg`, `rfortran-arpack`, BLAS, LAPACK, or `fortran-lapack` is vendored in this package. The current slice has no external Fortran dependency.
