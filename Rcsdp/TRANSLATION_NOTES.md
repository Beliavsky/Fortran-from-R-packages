# Translation notes

## Source

- R package: Rcsdp 0.1.57.6
- Bundled solver: CSDP 6.1.1
- Original license: Common Public License Version 1.0

## Module map

| Original code | Modern Fortran translation |
|---|---|
| `blockmat.h`, R conversion structures | `rcsdp_types.f90`, `rcsdp_problem.f90` |
| `Fnorm.c`, `norms.c`, `trace_prod.c` | `rcsdp_block_ops.f90`, `rcsdp_linalg.f90` |
| `allocmat.c`, `copy_mat.c`, `zero_mat.c`, `make_i.c`, `add_mat.c`, `addscaledmat.c`, `sym_mat.c` | allocatable assignment plus `rcsdp_block_ops.f90` |
| `mat_mult.c` | `mat_mult` in `rcsdp_block_ops.f90` |
| `chol.c`, `solvesys.c` | `rcsdp_linalg.f90`, solver Schur factorization |
| `op_a.c`, `op_at.c` | reference path in `rcsdp_problem.f90`; sparse path in `rcsdp_sparse_ops.f90` |
| `op_o.c` | reference path in `rcsdp_problem.f90`; sparse/hybrid path in `rcsdp_sparse_ops.f90` |
| `makefill.c` | `build_fill_workspace` in `rcsdp_fill_ops.f90` |
| `mat_multsp.c` | `mat_multspa`, `mat_multspb`, `mat_multspc` in `rcsdp_fill_ops.f90` |
| `easysdp.c` block classification / `nextbyblock` | `csdp_sparse_workspace` in `rcsdp_sparse_ops.f90` |
| `linesearch.c` | exact/Lanczos hybrid in `line_search_pd` |
| `qreig.c` | small LAPACK eigensolve of the Lanczos tridiagonal projection |
| Schur scaling/retry in `sdp.c` | `factor_schur` in `rcsdp_solver.f90` |
| Schur iterative refinement in `sdp.c` | `refine_schur_local` in `rcsdp_solver.f90` |
| `sdp.c`, `easysdp.c` | `csdp` in `rcsdp_solver.f90` |
| `readprob.c`, `writeprob.c`, `readsol.c`, `writesol.c` | `rcsdp_io.f90` |
| R `sparse.R` | `rcsdp_triplet.f90` |
| R `.Call`, S3/S4, `Matrix` coercions | intentionally omitted |

## v0.2 sparse Schur layer

v0.2 translated the key sparse Schur assembly strategy from CSDP. Constraint blocks are cross-indexed by SDP block, classified using CSDP's original sparse/dense rule, and contracted directly without permanently materializing all `A_i` matrices.

The Schur operator remains

`O(d) = A(Z^{-1} A'(d) X)`.

Sparse/sparse PSD pairs use the same four symmetry cases as CSDP `op_o.c`; diagonal blocks are matched directly by index. Dense-classified constraint blocks are materialized only block-by-block and use dense matrix products.

## v0.3 fill-restricted products

CSDP's `makefill` records every matrix entry that can be relevant to later calls to `A(.)`: all diagonal positions, nonzeros of `C`, and symmetric nonzero positions from every constraint.

`csdp_fill_workspace` stores that set per block and additionally builds:

- `row_ptr`, for products where the left operand is fill-restricted (`spA`);
- `col_ptr` plus `col_entry`, for products where the right operand is fill-restricted (`spB`);
- the explicit `(i,j)` fill list, for products where only selected output entries are required (`spC`).

The three translated kernels implement the same computational contracts as CSDP:

- `mat_multspa`: exploit fill in `A`;
- `mat_multspb`: exploit fill in `B`;
- `mat_multspc`: compute only output entries in fill.

For a block whose fill density exceeds `fill_density_limit` (default 0.01), the implementation falls back automatically to full `matmul`, matching the intent of CSDP's `SPARSELIMA/B/C` thresholds.

The predictor now preserves the `Zi*Fd` result from affine RHS formation and constructs

`dX = -(I - Zi*Fd + Zi*A'(dy))*X`,

which follows the original CSDP dataflow and removes redundant full matrix products. Corrector products and true-operator Schur refinement similarly use `spA`/`spC` where only fill entries can affect the next sparse operator.

## OpenMP

The fill data are grouped so parallel iterations write disjoint rows, columns, or individual output entries. `!$omp` directives are therefore included on the safe restricted-product loops. They are comments in a normal serial build and become active with compiler OpenMP flags such as `-fopenmp`.

Compiler-specific C prefetch builtins were not carried over; the Fortran implementation relies on contiguous row/column traversal and the compiler/runtime.

## Schur scaling and regularization

Before Cholesky, v0.3 preserves the unscaled Schur diagonal, computes its Euclidean norm, adds a CSDP-style retry diagonal when necessary, and forms the diagonal equilibration

`D(i) = 1/sqrt(O(i,i))`.

The factorized system is `D O D`; right-hand sides and solutions are correspondingly scaled by `D`, matching the role of CSDP's `workvec8`.

`use_schur_scaling=.false.` disables equilibration while retaining robust diagonal retry behavior for regression comparison.

## Large-block Lanczos line search

For a matrix block, the step boundary is controlled by the maximum eigenvalue of

`-R^{-T} dX R^{-1}`,

where `R` is the upper Cholesky factor of the current positive-definite block.

For blocks up to `lanczos_threshold` (default 180), v0.3 continues to compute all eigenvalues with LAPACK. Larger blocks use a fully reorthogonalized Lanczos iteration, capped at `lanczos_iterations` (default 30). The projected tridiagonal is only a few dozen rows, so its eigenvalues are obtained with the existing LAPACK symmetric eigensolver rather than reproducing CSDP's historical `qreig.c` routine.

A small conservative inflation is applied to the Ritz maximum because finite Lanczos iteration can underestimate the true limiting eigenvalue. The solver's later positive-definiteness safeguards remain in place.

## Compatibility paths

The performance layers can be disabled independently:

- `use_fill_products=.false.`: v0.2-style predictor/corrector products;
- `use_sparse_schur=.false.`: v0.1 dense-reference constraint and Schur path;
- `use_schur_scaling=.false.`: unscaled Schur Cholesky;
- `use_lanczos_linesearch=.false.`: exact LAPACK line search for all blocks.

These switches are intentionally retained for numerical regression and profiling.

## Validation

The new restricted products are compared directly with complete dense matrix multiplication. `spC` is checked only on the entries its contract promises to generate. `apply_o_fill` is checked against the v0.2 true sparse operator.

The complete solver is also run in default v0.3 mode, no-fill mode, and no-Schur-scaling mode on the bundled `theta1` problem. All converge to the same objective within the regression tolerance. A separate 120x120 SDP with diagonal-only fill verifies that the integrated predictor/corrector path actually exercises the sparse restricted kernels rather than their dense fallback.

The Lanczos line search is tested on a 220x220 problem with a known limiting eigenvalue and compared with the exact LAPACK path.

## Remaining low-level differences

- packed C work arrays are replaced by allocatable Fortran arrays;
- C compiler prefetch intrinsics are omitted;
- OpenMP scheduling is portable and simpler than every original C pragma;
- Lanczos uses LAPACK on the projected tridiagonal rather than Algorithm 464 `qreig.c`;
- some historical stopping and refinement heuristics are expressed more directly;
- the direct Schur matrix remains dense, consistent with CSDP's architecture.

## License

The sparse Schur, fill, restricted-product, scaling, and line-search algorithms translated here derive from CSDP 6.1.1 sources distributed under CPL-1.0. The original license and provenance files remain in the project unchanged.
