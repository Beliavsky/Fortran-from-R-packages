# sanic — modern Fortran translation

This FPM package translates the computational API of the R package
[`sanic`](https://github.com/nk027/sanic), version 0.0.2, to modern free-form
Fortran. It provides direct linear-system solvers, sparse Krylov solvers,
eigendecomposition, SVD, and Arnoldi/Lanczos routines without Rcpp, Eigen,
BLAS/LAPACK system links, or other separately installed numerical libraries.

The translation preserves the upstream GPL (>= 3) notice and provenance. See
`LICENSE`, `NOTICE.md`, `COMPATIBILITY.md`, and `upstream/`. The preserved
Spectra header is provenance-only and remains under MPL-2.0.

## Build

Extract this project as a top-level directory named `sanic` beside any other
repositories you keep under `Fortran-from-R-packages`, then run:

```text
fpm build
fpm test
fpm run --example basic_solve
```

No system `-lblas` or `-llapack` link is used.

The maintained Fortran sources use one real kind, `dp = real64`, exported by
the public `sanic` module. All newly translated code is free-form Fortran.

## Public API

The umbrella module is `sanic`. Important procedures and types include:

- `sparse_matrix`, `sparsify`, `sparse_from_triplets`, `sparse_to_dense`
- `solve2_dense`, `solve2_sparse`
- `solve_chol_dense`, `solve_chol_sparse`
- `solve_lu_dense`, `solve_lu_sparse`
- `solve_qr_dense`, `solve_qr_sparse`
- `solve_cg_sparse` with `sanic_bicgstab`, `sanic_lscg`, and `sanic_cg`
- `eigen2` and `eigen_result`
- `svd2`, `svd_result`, `sanic_svd_bdc`, and `sanic_svd_jacobi`
- `arnoldi_sparse`, `lanczos_sparse`, and `krylov_result`

R-like generic names `solve2`, `solve_chol`, `solve_lu`, and `solve_qr` are also
exported from the umbrella module for dense/sparse overload dispatch.

### Example

```fortran
program basic_solve
   use sanic, only: dp, solve2_dense
   implicit none

   real(dp) :: a(2, 2), b(2, 1)
   real(dp), allocatable :: x(:, :)
   integer :: info

   a = reshape([4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
   b(:, 1) = [6.0_dp, 7.0_dp]
   call solve2_dense(a, x, info, b=b)
   if (info /= 0) error stop "solve failed"
   write(*, '(a,2(1x,f10.6))') "x =", x(:, 1)
end program basic_solve
```

Expected result: `x = 1, 2`.

## Numerical behavior

Dense direct paths use package-local LLT/pivoted LDLT, Gaussian LU, and
Householder QR algorithms. Dense QR now keeps the three upstream selector
behaviors distinct: unpivoted Householder QR, column-pivoted Householder QR,
and complete row/column pivoting. Underdetermined solves return a valid
rank-revealing basic solution, matching Eigen QR solve semantics more closely
than the earlier pseudoinverse/minimum-norm fallback.

Sparse Cholesky and LU operate on dynamic sparse rows and permit fill-in without
first expanding the coefficient matrix to dense storage. Natural ordering is
implemented directly. The `ordering=0` and `ordering=1` paths use dynamic
minimum-degree elimination on sparse adjacency lists built from symmetrized or
column-intersection graphs. These remain structural approximations to
AMD/COLAMD rather than exact SuiteSparse/Eigen orderings. Sparse QR remains
sparse during numerical factorization, using Givens rotations, structural
ordering, and numerical column pivoting for tall, square, and wide systems.
Sparse symmetry detection canonicalizes duplicate coordinates and uses binary
transpose lookup rather than dense storage or repeated full triplet scans.

Sparse CG, BiCGSTAB, LSCG, Arnoldi, and Lanczos use sparse matrix-vector
products directly. `precond=0` selects identity preconditioning and
`precond=1` Jacobi/diagonal preconditioning. For `precond=2`, BiCGSTAB uses a
package-local dual-threshold ILUT with fill limiting and fill-reducing
permutation. CG uses a scaled, permuted, adaptively shifted limited-memory
incomplete Cholesky that retains the largest generated fill entries up to the
original column capacity. LSCG intentionally keeps the upstream/default
diagonal preconditioner. The default iteration cap follows Eigen's solver
default of `2*ncol` when the caller does not supply a positive value. Exact
bitwise or iteration-for-iteration identity with Eigen is not claimed.

Symmetric eigenproblems primarily use Householder tridiagonalization followed
by implicit QL iteration, with the earlier Jacobi eigensolver retained as a
fallback. When `symmetric=.true.` is supplied, only the diagonal and lower
triangle define the self-adjoint matrix, matching Eigen's `SelfAdjointEigenSolver`.
General real eigenproblems first reduce to upper Hessenberg form,
then use shifted complex Householder-QR iteration and inverse iteration for
right eigenvectors. Complex eigenvalues are supported.

The two SVD selectors now use distinct numerical paths. `BDC` uses Householder
bidiagonalization followed by implicit bidiagonal QR iteration; this reproduces
the bidiagonal numerical structure but is not Eigen's divide-and-conquer BDC
implementation. `Jacobi` uses column-pivoted Householder QR preconditioning and
two-sided Jacobi sweeps on the square core, with the direct one-sided Jacobi
routine retained as a convergence fallback. Both paths support values-only,
thin-vector, and full-vector output for tall and wide matrices without forming
a Gram matrix.

For Arnoldi and Lanczos, omitting the starting vector uses the deterministic
vector `1,2,...,n`, normalized internally. Upstream R instead calls `rnorm()`.
On early Krylov breakdown, result arrays preserve the requested `iter` shapes
and unused entries are deterministically zero-filled; `iterations` records the
achieved Krylov dimension.

### Literal upstream defects intentionally not reproduced

A source-level compatibility audit found several defects in sanic 0.0.2 that
would make literal bug-for-bug compatibility conflict with the package's own
documented API. Solver objects declared inside C++ selector branches shadow the
objects used after those branches, so several pivot, ordering, and
preconditioner choices are accidentally ignored upstream. The R `solve_cg()`
wrapper also validates/coerces an explicit `iter` from `tol` rather than from
`iter`, and its `length(b)` dimension check rejects multi-column RHS matrices.
The Fortran translation implements the documented selector behavior, honors
`iter`, and supports multiple dense RHS columns. See `COMPATIBILITY.md` for the
exact source locations and rationale.

## Translation coverage

Package status: **substantial**.

**10 of 10 (100%)** exported computational R functions have meaningful Fortran
mappings. The fraction measures mapped computational functions, not complete R
interface compatibility or numerical identity with Eigen.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `sparsify` | `sanic_sparse` | `sparsify` | substantial |
| `solve2` | `sanic_direct` | `solve2_dense`, `solve2_sparse` | substantial |
| `solve_chol` | `sanic_direct` | `solve_chol_dense`, `solve_chol_sparse` | substantial |
| `solve_lu` | `sanic_direct` | `solve_lu_dense`, `solve_lu_sparse` | substantial |
| `solve_qr` | `sanic_direct` | `solve_qr_dense`, `solve_qr_sparse` | substantial |
| `solve_cg` | `sanic_iterative` | `solve_cg_sparse` | substantial |
| `eigen2` | `sanic_spectral` | `eigen2` | partial |
| `svd2` | `sanic_spectral` | `svd2` | partial |
| `arnoldi` | `sanic_krylov` | `arnoldi_sparse`, `lanczos_sparse` | substantial |
| `lanczos` | `sanic_krylov` | `lanczos_sparse` | substantial |

See `API_COVERAGE.md` and the `[extra.translation]` entries in `fpm.toml` for
per-function notes.

## Tests

Deterministic FPM tests cover:

- partial/full-pivot LU and inverse calculation;
- LLT and LDLT symmetric solves, including sparse fill-in cases;
- Householder QR for full-rank, rank-deficient, tall, wide basic-solution, and sparse-ordered systems;
- sparse storage, duplicate triplets, matrix-vector/matrix products, and sparse symmetry checks;
- CG, BiCGSTAB, and least-squares CG, including dual-threshold ILUT, limited-memory incomplete Cholesky, and LSCG `precond=2` behavior;
- symmetric real and nonsymmetric complex eigenproblems, including strongly nonnormal eigenvector residuals;
- distinct bidiagonal-QR and Jacobi SVD paths, thin/full/values-only shapes, wide matrices, reconstruction, and singular values down to `1e-12`;
- Lanczos and Arnoldi projected eigenvalues plus requested-shape preservation on early breakdown;
- Eigen-compatible lower-triangle semantics for explicitly self-adjoint eigensolves;
- invalid selector/tolerance rejection, explicit iterative budgets, and multiple iterative RHS columns.

The tests use no random-number generator.
