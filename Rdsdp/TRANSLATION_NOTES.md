# Translation notes

## Scope

This project translates the computational functionality used by **Rdsdp
1.0.6** from its bundled DSDP5 source tree into modern Fortran.  The embedded
C build metadata contains paths naming DSDP5.8.  The translation is
algorithmic rather than ABI-compatible: native Fortran derived types replace
R `.Call` objects and C pointer layouts.

v0.3.0 extends the v0.2 sparse/low-rank data layer with two Schur-system
algorithms: reusable sparse direct factorization and matrix-free CG.

## Source mapping

| Original area | Fortran implementation |
|---|---|
| `R/Rdsdp.R`, `src/Rdsdp.c` semantics | `src/rdsdp.f90`, `src/rdsdp_problem.f90` |
| DSDP solver driver / dual algorithm | `src/rdsdp_solver.f90` |
| `solver/dsdpcg.c` role | assembled and matrix-free PCG in `src/rdsdp_solver.f90` |
| `sdp/dsdpdatamat*.c`, `sdp/dsdpadddatamat.c` | `src/rdsdp_data.f90`, `src/rdsdp_types.f90` |
| SDP cone equations | SDP branches in `src/rdsdp_solver.f90` |
| LP cone equations | LP branches in `src/rdsdp_solver.f90` |
| DSDP sparse Schur direct-solve role (`sdpsymb`, `sdporder`, `sdpnfac`, `cholmat*`) | `src/rdsdp_sparse_ldlt.f90` using RCM + QDLDL-derived sparse LDL^T |
| dense matrix/Cholesky work | `src/rdsdp_linalg.f90` + BLAS/LAPACK |
| SDPA reader | `src/rdsdp_io.f90` |
| option-file semantics | `src/rdsdp_options.f90` |

The v0.3 sparse direct module provides the same category of capability as the
historical DSDP custom sparse factorization stack, but it is not presented as
a line-by-line translation of those particular files.

## Mathematical conventions retained

For an SDP block:

```text
S_k = C_k - sum_i y_i A_{k,i} + r I
```

and for a linear/nonnegative block:

```text
s = c - A y + r 1.
```

The residual shift `r > 0` implements the infeasible start.  Primal recovery
retains the DSDP Newton correction

```text
X = mu*S^{-1} + mu*S^{-1} A^*(d) S^{-1}.
```

## Sparse/low-rank data layer from v0.2

The SDP data abstraction supports:

```text
dsdp_data_dense
dsdp_data_sparse
dsdp_data_lowrank
```

Sparse matrices use full-symmetric COO.  Low-rank matrices use weighted outer
products.  SDPA sparse input is allocated directly in sparse form.

The SDP Schur term is

```text
H_ij = mu * trace(S^{-1} A_i S^{-1} A_j).
```

Sparse-sparse, sparse-low-rank, and low-rank-low-rank contractions avoid an
`n x n x m` transformed tensor.

## v0.3 matrix-free Hessian action

For an SDP block and a Newton-space vector `p`, define

```text
B = sum_i p_i A_i - p_r I
T = S^{-1} B S^{-1}.
```

Then the block contribution to the Hessian-vector product is

```text
(Hp)_i += mu * <A_i,T>
(Hp)_r -= mu * trace(T).
```

The residual barrier separately contributes `mu*p_r/r^2`.

For an LP block with `w2 = s^{-2}`:

```text
u = A p_y - p_r 1
(Hp)_y += mu * A^T (w2 .* u)
(Hp)_r -= mu * sum(w2 .* u).
```

The exact diagonal is built during slack evaluation and used as a Jacobi
preconditioner.  This removes the `m x m` Schur allocation/assembly from CG
iterations.  A failed matrix-free solve can assemble the normal Schur matrix
and use the existing fallback sequence.

## v0.3 sparse direct Schur solve

`src/rdsdp_sparse_ldlt.f90` implements:

- structural/numerical sparsity extraction from the assembled symmetric Schur
  matrix;
- reverse Cuthill-McKee ordering;
- upper-triangular CSC storage;
- elimination-tree symbolic analysis;
- QDLDL-derived sparse LDL^T numerical factorization;
- forward/diagonal/back substitution;
- cached symbolic pattern reuse across Newton iterations.

The QDLDL-derived kernel is licensed under Apache-2.0.  The RCM ordering and
cache/reuse integration are part of this translation.

With default `sparse_schur_drop_tol=0`, no numerically small nonzero is dropped;
only exact zeros are omitted.  This keeps the sparse path conservative.

The sparse direct path is attempted only when the Schur dimension and density
meet the user controls.  Dense or failed sparse systems fall back safely.

## Why both v0.2 and v0.3 paths remain

Sparse data do not imply a sparse Schur matrix.  In a single SDP block,
`S^{-1}` is typically dense, so most constraint pairs can interact.  The
`mcp100` case is a good example: direct sparse pair contractions make Schur
assembly cheap, while the resulting Schur system is not sparse enough to
benefit from the v0.3 sparse direct factor.

Conversely, LPs or multi-block SDPs with weakly coupled constraint groups can
produce genuinely sparse Schur matrices, where symbolic reuse and sparse
factorization are useful.

Matrix-free CG addresses a different regime: large `m` where the Schur matrix
may be dense but Hessian-vector products are much cheaper than building and
factoring all pairwise entries.

## Validation specific to v0.3

The v0.3 suite adds:

- direct sparse LDL^T residual tests;
- a changed-diagonal refactorization test proving one symbolic analysis is
  reused for two numerical factorizations;
- a 120-constraint LP whose arrowhead Schur systems use the sparse factor path
  throughout with no dense fallback;
- a dense-constraint orthogonal LP that forces matrix-free PCG and converges
  without a direct fallback.

All earlier R examples, SDPA tests, original DSDP C objective references,
low-rank tests, sparse-data tests, and `mcp100` remain active.

## Numerical linear algebra

BLAS/LAPACK remain required for SDP slack Cholesky/inversion, line-search
log-determinants, primal recovery, and the universal dense Schur fallback.
Sparse LDL^T is self-contained Fortran except for the rest of the solver's
normal BLAS/LAPACK dependencies.

## Remaining differences

- The sparse direct implementation is RCM + QDLDL-derived LDL^T, not an exact
  reproduction of DSDP's custom sparse Cholesky source files.
- Matrix-free PCG currently uses Jacobi preconditioning; richer historical
  DSDP preconditioner/reuse modes are not all reproduced.
- Automatic low-rank detection, specialized update/eigen kernels, extended
  bound-cone APIs, and platform-specific parallel/cache tuning remain omitted.
- Full certificate behavior for every historical infeasible/unbounded exit is
  not claimed.

## License preservation

- Rdsdp GPL-3 material is retained under `licenses/`.
- The DSDP University of Chicago/Argonne notice is retained in `LICENSE` and
  `licenses/DSDP-LICENSE`.
- Apache-2.0 for the QDLDL-derived sparse factor kernel is retained in
  `licenses/QDLDL-APACHE-2.0.txt`.
