# Changelog

## 0.3.0

- Added `rcsdp_fill_ops`, translating CSDP `makefill.c` and `mat_multsp.c`.
- Added reusable row/column fill indices and the `mat_multspa`, `mat_multspb`, and `mat_multspc` kernels.
- Integrated fill-restricted products into affine RHS formation, predictor reuse, corrector RHS/direction formation, and true-operator Schur refinement.
- Added portable `!$omp` parallel loops for race-free restricted-product traversals; serial builds remain unchanged.
- Added CSDP-style Schur diagonal equilibration and diagonal retry regularization.
- Added a large-block Lanczos line-search path with double reorthogonalization and conservative Ritz scaling; exact LAPACK eigenvalue search remains the small-block/default fallback.
- Added controls `use_fill_products`, `fill_density_limit`, `use_schur_scaling`, `use_lanczos_linesearch`, `lanczos_threshold`, and `lanczos_iterations`.
- Added solution diagnostics for fill size, sparse/dense fill-product counts, Lanczos use, and maximum Schur diagonal regularization.
- Added `test_fill_products`, `test_solver_modes`, `test_sparse_fill_solver`, and `test_lanczos_linesearch`.
- Added `benchmark_fill_products`; a local 500x500/0.296%-fill run showed about a 44x speed ratio for `spC` versus a full dense product.
- Preserved the v0.2 no-fill path and v0.1 dense Schur path for A/B regression testing.

## 0.2.0

- Added `rcsdp_sparse_ops`, a native Fortran translation of the important sparse constraint/Schur kernels from CSDP 6.1.1.
- Added block-wise cross-indexing equivalent to CSDP's `nextbyblock`/`byblocks` structure.
- Added CSDP's original sparse/dense constraint-block classification rule.
- Added direct sparse-entry `A(X)` and `A'(y)` operators.
- Replaced the default dense-column Schur construction with sparse-aware pairwise assembly.
- Added hybrid dense-block Schur handling for constraint blocks classified as dense.
- Avoided persistent dense copies of all constraint matrices in the default solver path.
- Added true-operator Schur residual checking and iterative refinement when `fastmode=.false.`.
- Added `use_sparse_schur` control; `.false.` retains the v0.1 dense-reference implementation.
- Added solution statistics for constraint nonzeros, sparse/dense block counts, Schur assemblies, pair contractions, dense products, refinements, and Schur CPU time.
- Added `test_sparse_schur`, comparing sparse and dense `A`, `A'`, and Schur operators on the bundled `theta1` problem.
- Added dense-fallback solver regression coverage.
- Added `benchmark_theta` example comparing the v0.2 sparse path with the v0.1 dense-reference path.
- Preserved the existing public SDP problem/solution API and CPL-1.0 licensing/provenance.

## 0.1.0

- Initial modern Fortran/FPM translation of Rcsdp 0.1.57.6 / CSDP 6.1.1 computational code.
- Added block-diagonal SDP/LP data structures and sparse symmetric constraint storage.
- Ported CSDP predictor-corrector Newton equations, initialization, stopping measures, controls, and status codes.
- Added BLAS/LAPACK-backed Cholesky, triangular inversion, Schur solves, and exact eigenvalue line search.
- Added SDPA sparse problem and solution I/O.
- Added Rcsdp-style symmetric triplet matrix helpers.
- Added regression tests for the canonical CSDP example, diagonal LP, manifold example, I/O round trips, and bundled theta1 problem.
- Preserved CPL-1.0 license and original authorship/provenance files.
