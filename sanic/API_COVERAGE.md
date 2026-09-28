# API coverage

Coverage basis: distinct exported computational R functions in `NAMESPACE`.
The upstream package exports 10 such functions, and all 10 have a meaningful
Fortran mapping. This is **function mapping coverage**, not a claim of complete
R, Eigen, or sparse-backend compatibility.

| R function | Fortran module | Public Fortran procedure(s) | Status |
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

## Status rationale

Overall status is **substantial**. The complete exported computational surface
has been mapped. Dense direct paths include pivoted LU, LLT/pivoted LDLT, and
three distinct Householder QR selector paths, including rank-revealing basic
solutions for underdetermined systems. Sparse Cholesky, LU, and QR operate on
dynamic sparse rows without densifying the coefficient matrix; sparse symmetry
checking canonicalizes duplicate coordinates and performs binary transpose
lookups. Sparse iterative methods operate directly on sparse storage. The
upstream `precond=2` BiCGSTAB path uses dual-threshold ILUT with fill limiting,
while CG uses scaled, permuted, shifted limited-memory incomplete Cholesky.

The spectral paths are also distinct: symmetric eigenproblems use
tridiagonalization plus implicit QL with a Jacobi fallback; general
eigenproblems use Hessenberg reduction plus shifted complex QR; explicit
self-adjoint solves use the lower triangle as Eigen does; BDC-selected SVD uses
bidiagonalization plus implicit QR, while Jacobi-selected SVD uses
column-pivoted QR preconditioning and two-sided Jacobi sweeps. Arnoldi/Lanczos
preserve requested output shapes on early breakdown while separately reporting
the achieved Krylov dimension.

The translation intentionally follows sanic's documented selector behavior
rather than reproducing C++ block-scope shadowing defects that cause several
upstream pivot, ordering, and preconditioner branches to be ignored. It also
honors the caller's explicit iterative `iter` value rather than the upstream R
wrapper typo that coerces `tol` into `iter`. These differences are documented
in `COMPATIBILITY.md`.

The package is not marked complete because its dynamic minimum-degree
orderings approximate rather than reproduce Eigen/SuiteSparse AMD/COLAMD
exactly, incomplete-factorization drop/fill details are not guaranteed to be
numerically identical to Eigen, the `BDC` path is bidiagonal QR rather than
Eigen's divide-and-conquer implementation, and the portable eigensolvers need
not reproduce Eigen's eigenvalue ordering or roundoff behavior exactly.

The exact machine-readable mapping and compatibility notes live in
`fpm.toml` under `[extra.translation]` and `[[extra.translation.function]]`.
