# Rdsdp-fortran

`Rdsdp-fortran` is a modern Fortran/FPM translation of the computational path
exposed by **Rdsdp 1.0.6** and its bundled **DSDP5** solver sources.  The
embedded source/build metadata identifies the bundled solver tree as DSDP5.8.

v0.3.0 keeps the v0.2 public API and adds the linear-system layer that was the
main remaining scalability boundary: reusable sparse Schur factorization for
structurally sparse Newton systems and matrix-free preconditioned CG for large
systems where forming the Schur matrix is undesirable.

## Implemented computational functionality

- Semidefinite (PSD) blocks and linear/nonnegative blocks.
- SeDuMi-style input compatible with `Rdsdp::dsdp()` through the Fortran
  `dsdp(A,b,C,l,s,...)` procedure.
- DSDP dual slack convention `S = C - A^*(y) + r I` with infeasible-start
  residual shift `r`.
- Dual logarithmic-barrier / dual-scaling Newton iterations.
- DSDP-style primal recovery, including the Newton correction to
  `mu*S^{-1}` required for accurate primal feasibility.
- Dense, sparse symmetric COO, and weighted low-rank SDP data matrices.
- Sparse-aware Schur assembly using direct sparse/low-rank contractions.
- Direct LAPACK Cholesky Schur solve as the universal fallback.
- **RCM-ordered sparse LDL^T Schur solve** with symbolic-pattern reuse.
- **Matrix-free Jacobi-preconditioned CG**, evaluating Hessian-vector products
  from the current cone slacks without materializing the full Schur matrix.
- Assembled-matrix CG retained for comparison with v0.2 behavior.
- Sparse SDPA-format input/output and DSDP-style option-file parsing.
- Solver statistics exposing which data, factorization, and iterative paths
  were actually used.

## Sparse and low-rank SDP data

For SDP blocks the supported storage kinds are:

- `dsdp_data_dense`
- `dsdp_data_sparse`
- `dsdp_data_lowrank`

Sparse SDPA input is read directly into full-symmetric COO storage.  Dense
SeDuMi-style input is compressed matrix-by-matrix when its density is below
`control%sparse_density_threshold` (default 0.20).

Low-rank matrices use

```text
A = sum_j coeff(j) v_j v_j^T
```

and can be installed with `set_constraint_lowrank` and
`set_objective_lowrank`.

## v0.3 sparse Schur factorization

For assembled Newton matrices that are genuinely sparse, v0.3.0 can avoid the
dense LAPACK factorization.  The path:

1. identifies the numerical/structural Schur pattern;
2. computes a reverse Cuthill-McKee ordering;
3. stores the permuted upper triangle in CSC form;
4. performs elimination-tree symbolic analysis;
5. factors with a native sparse LDL^T algorithm derived from QDLDL; and
6. reuses the symbolic analysis while only the numerical values change.

The controls are:

```fortran
ctrl%use_sparse_schur_factor = .true.   ! default
ctrl%sparse_schur_threshold = 80
ctrl%sparse_schur_density_limit = 0.20_dp
ctrl%sparse_schur_drop_tol = 0.0_dp
```

A zero drop tolerance is the conservative default: only exact structural zeros
are omitted.  If the matrix is too dense or the sparse factorization fails,
the solver falls back to the other enabled Schur path and ultimately LAPACK
Cholesky.

This is functionally analogous to the sparse direct role of DSDP's historical
`sdpsymb/sdporder/sdpnfac/cholmat*` stack, but it is **not** a line-for-line
translation of those files.  The sparse numerical kernel is QDLDL-derived and
therefore carries Apache-2.0 licensing in addition to the existing Rdsdp/DSDP
licenses.

## v0.3 matrix-free CG

When CG is enabled, matrix-free mode is now the default CG implementation:

```fortran
ctrl%use_cg = .true.
ctrl%cg_matrix_free = .true.
ctrl%cg_threshold = 200
ctrl%cg_tol = 1.0e-10_dp
ctrl%cg_maxiter = 80
ctrl%cg_fallback_direct = .true.
```

For an SDP block, a Hessian-vector product forms

```text
B = sum_i p_i A_i - p_r I
T = S^{-1} B S^{-1}
```

and contracts `T` with each constraint matrix.  For an LP block it evaluates
`A^T diag(s^{-2}) A p` through two matrix-vector products.  Thus no `m x m`
Schur matrix is required for the CG iterations.

A Jacobi preconditioner is built from the exact Schur diagonal.  If matrix-free
CG does not meet its residual target and `cg_fallback_direct` is true, the
assembled Schur matrix is formed and the normal fallback sequence is used.

Matrix-free CG is not universally faster.  On problems such as `mcp100`, where
individual SDP constraint matrices are extraordinarily sparse, the v0.2 direct
pair-contraction assembly remains faster.  Matrix-free mode is most useful
when `m` is large and forming/factoring a dense Schur matrix dominates.

Set

```fortran
ctrl%cg_matrix_free = .false.
```

to retain the assembled-matrix PCG path from v0.2.

## Public API

The top-level module is `rdsdp`.

```fortran
use rdsdp

real(dp), allocatable :: a(:,:), b(:), c(:)
integer :: l
integer, allocatable :: s(:)
type(dsdp_solution) :: sol
type(dsdp_control) :: ctrl

! ... fill a, b, c, l, s ...
ctrl = dsdp_control()
call dsdp(a, b, c, l, s, sol, ctrl)
print *, sol%dobj, sol%y
```

For an SDPA sparse file:

```fortran
use rdsdp

type(dsdp_solution) :: sol
call dsdp_readsdpa('problem.dat-s', sol)
print *, sol%dobj, sol%pobj
```

The lower-level `read_sdpa`, `write_sdpa`, `dsdp_solve`,
`dsdp_from_sedumi`, `get_data_dense`, low-rank setters, and the
problem/control/solution derived types remain public.

## Diagnostics

In addition to the v0.2 counters, `dsdp_solution` now reports:

- `matrix_free_cg_solves`
- `matrix_free_matvecs`
- `sparse_factor_solves`
- `sparse_factor_fallbacks`
- `sparse_symbolic_analyses`
- `sparse_numeric_factorizations`
- `schur_matrix_nnz`
- `schur_factor_nnz`

These make symbolic reuse and fallback behavior directly testable.

## Option-file additions

The option parser accepts the v0.3 keys:

```text
-cgmatrixfree
-usesparsefactor
-sparsefactorthreshold
-sparsefactordensity
-sparsefactordroptol
```

along with the v0.1/v0.2 DSDP-style controls.

## Building with FPM

BLAS and LAPACK are required.

```text
fpm build
fpm test
fpm run --example solve_sdpa -- data/control1.dat-s
fpm run --example benchmark_mcp100
fpm run --example benchmark_matrixfree_lp
fpm run --example benchmark_sparse_factor
```

The supplied `fpm.toml` links `lapack` and `blas`.  FPM was not installed in
the translation environment, so the release was validated directly with GNU
Fortran 14.2.0 using the same source/module dependency order.

## Validation

The original R-package regression solutions continue to reproduce:

- approximately `[-1.0, -0.75]`
- approximately `[1.0]`
- approximately `[0.6, -0.4, 3.0]`

The independent bundled-DSDP C reference checks remain within their established
tolerances:

| problem | bundled DSDP C reference | Fortran v0.3 error |
|---|---:|---:|
| `control1.dat-s` | -17.7846267761 | about 1.0e-6 |
| `truss1.dat-s` | 8.99999631296 | about 7.1e-7 |
| `vibra1.dat-s` | -40.8190124025 | about 2.4e-6 |

`mcp100.dat-s` continues to converge to approximately
`-226.15735223`.

The v0.3 tests additionally verify:

- sparse LDL^T residual accuracy;
- symbolic analysis reuse across changed numerical diagonals;
- an integrated 120-constraint arrowhead Schur solve with no dense fallback;
- matrix-free CG convergence without Schur materialization;
- all v0.2 sparse/low-rank/CG/dense-fallback tests;
- SDPA and option-file round trips.

## Local performance measurements

Measurements below use GNU Fortran 14.2.0, system BLAS/LAPACK, and `-O2` in
the translation environment.  They are illustrative rather than portable
benchmarks.

### Sparse direct factor benchmark

For a 1000 x 1000 tridiagonal SPD system:

| backend | CPU time | factor/storage indicator |
|---|---:|---:|
| RCM + sparse LDL^T | about 0.0168 s | 999 strict-lower factor nnz |
| dense LAPACK Cholesky | about 0.1982 s | 1,000,000 dense matrix entries |

The local dense/sparse time ratio was about **11.8x**, with residuals around
`5e-15` for both solves.

### Matrix-free CG benchmark

For the included generated 300-variable dense-constraint LP:

| backend | CPU time | objective |
|---|---:|---:|
| assembled dense Schur + direct solve | about 2.13 s | 299.999683 |
| matrix-free PCG | about 0.424 s | 299.999683 |

The local time ratio was about **5.0x**.  Matrix-free PCG required 63 Hessian
matvecs and no direct fallback.

### `mcp100`

The existing sparse-data benchmark remains approximately:

| backend | CPU time | objective |
|---|---:|---:|
| sparse data/Schur contraction | about 0.469 s | -226.15735223 |
| retained dense data path | about 23.72 s | -226.15735223 |

This demonstrates why v0.3 retains rather than replaces the v0.2 pairwise
sparse Schur assembly path.

## Remaining differences from original DSDP5

Major exposed Rdsdp computations are now present, but several implementation
and extended-API differences remain:

- the sparse direct solver uses RCM + QDLDL-derived LDL^T rather than DSDP's
  historical custom ordering/Cholesky implementation;
- matrix-free CG uses a Jacobi preconditioner rather than every original DSDP
  factored/unfactored preconditioner mode;
- automatic low-rank detection/factorization of arbitrary input matrices is
  not implemented;
- specialized low-rank update/eigen kernels and platform-specific cache,
  prefetch, and OpenMP tuning are omitted;
- extended bound-cone APIs beyond the ordinary Rdsdp `K$l` path are omitted;
- full infeasibility/unboundedness certificate machinery for every historical
  DSDP exit path is not reproduced.

## Deliberately omitted interface code

- R `.Call` registration/marshalling and S3/S4 integration.
- R `methods`/`utils` glue.
- C pointer ownership and allocator wrappers.
- Plotting code (there is no material plotting algorithm in the solver path).

## Licensing and provenance

Three licensing/provenance layers are retained:

1. **Rdsdp 1.0.6:** GPL-3, with DESCRIPTION/NAMESPACE and GPL text under
   `licenses/`.
2. **DSDP5:** University of Chicago/Argonne permission notice, retained in
   `LICENSE` and `licenses/DSDP-LICENSE`.
3. **QDLDL-derived sparse LDL^T kernel:** Apache License 2.0, retained as
   `licenses/QDLDL-APACHE-2.0.txt`.

See `TRANSLATION_NOTES.md` for source mapping and implementation details.
