# Changelog

## 0.3.0

Schur linear-system architecture release.

### Added

- `src/rdsdp_sparse_ldlt.f90` with reverse Cuthill-McKee ordering, CSC Schur
  storage, elimination-tree symbolic analysis, QDLDL-derived sparse LDL^T
  factorization, triangular solves, and reusable symbolic-factor cache.
- Automatic sparse-direct selection using Schur dimension and density limits,
  with dense fallback.
- Matrix-free Hessian-vector products for SDP and LP blocks.
- Matrix-free Jacobi-preconditioned CG that avoids allocating or assembling an
  `m x m` Schur matrix during iterative solves.
- Assembled PCG remains selectable with `cg_matrix_free=.false.`.
- New controls: `cg_matrix_free`, `use_sparse_schur_factor`,
  `sparse_schur_threshold`, `sparse_schur_density_limit`, and
  `sparse_schur_drop_tol`.
- Option-file keys for all v0.3 controls.
- New diagnostics for matrix-free matvecs/solves, sparse-factor solves and
  fallbacks, symbolic/numeric factorization counts, and Schur/factor nnz.
- `test_sparse_factor` covering sparse residual accuracy, symbolic reuse, and
  a live arrowhead-Schur solver case.
- `test_matrixfree_cg` forcing the matrix-free path with no direct fallback.
- `benchmark_sparse_factor` and `benchmark_matrixfree_lp` examples.
- Apache-2.0 license text for the QDLDL-derived factorization kernel.

### Performance

Local GNU Fortran 14.2.0 `-O2` measurements:

- 1000 x 1000 tridiagonal SPD solve: about 0.0168 s with RCM+sparse LDL^T
  versus 0.1982 s with dense LAPACK Cholesky, about an 11.8x ratio.
- Generated 300-variable dense-constraint LP: about 0.424 s matrix-free PCG
  versus 2.13 s assembled-direct, about a 5.0x ratio, identical objective and
  no direct fallback.
- `mcp100`: sparse pair-assembly path remains preferable, about 0.469 s versus
  23.72 s for the retained dense-data reference path in the final local run.

### Compatibility

- v0.1/v0.2 public APIs and problem/data representations are unchanged.
- LAPACK Cholesky remains the universal fallback.
- The v0.2 sparse-pair Schur assembly path is retained and remains the default
  behavior on dense Schur systems such as `mcp100`.
- Sparse Schur factoring is enabled by default only above its threshold and
  below its density limit; CG remains opt-in.

### Remaining differences

- The sparse direct backend is RCM + QDLDL-derived LDL^T rather than an exact
  line-for-line port of DSDP's `sdpsymb/sdporder/sdpnfac/cholmat*` sources.
- Matrix-free CG uses Jacobi preconditioning rather than all historical DSDP
  factored/unfactored preconditioner modes.
- Automatic low-rank detection, extended bound cones, and platform-specific
  tuning remain deferred.

## 0.2.0

Sparse/low-rank data and Schur performance release.

### Added

- DSDP-style abstract SDP data storage with dense, symmetric sparse COO, and
  weighted low-rank outer-product representations.
- `src/rdsdp_data.f90` with matrix add/dot/materialization, sparse compression,
  low-rank setters, and sparse/low-rank Schur contractions.
- Two-pass SDPA reader that stores SDP matrices directly in sparse form rather
  than allocating a dense `n x n x m` constraint tensor.
- Automatic per-matrix compression of sparse SeDuMi-style input.
- Sparse-aware Schur assembly using direct data-matrix contractions.
- Optional diagonally preconditioned conjugate-gradient Schur solve with direct
  Cholesky fallback.
- Cholesky-only `spd_logdet` path for line-search states.
- Solver statistics for data nnz, Schur pair types, CG iterations, and Schur
  assembly/solve time.
- Public `set_constraint_lowrank` and `set_objective_lowrank` procedures.
- Sparse-backend, low-rank, CG, and dense-fallback regression tests.
- `mcp100.dat-s` full solve regression.
- `benchmark_mcp100` example.

## 0.1.0

Initial modern Fortran/FPM computational translation of Rdsdp 1.0.6 and its
bundled DSDP5 solver path.

### Added

- Native Fortran `dsdp(A,b,C,l,s,...)` interface.
- `dsdp_readsdpa` and sparse SDPA read/write support.
- SDP and linear/nonnegative cone blocks.
- DSDP-style infeasible-start residual shift.
- Dual logarithmic-barrier Newton/Schur iterations.
- Cholesky-based positive-definite line search.
- Corrected DSDP-style primal reconstruction.
- Objective, gap, feasibility, residual, and solution-type diagnostics.
- Common DSDP options-file parsing.
- FPM package metadata and BLAS/LAPACK linkage.
- Regression tests reproducing the three R-package examples.
- Independent objective regression against the bundled original DSDP C solver
  for `control1`, `truss1`, and `vibra1`.
- Preserved Rdsdp GPL-3 and DSDP license/provenance material.
