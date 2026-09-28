# Rlinsolve — modern Fortran translation

This directory translates the computational API of the R package **Rlinsolve
0.3.3** to modern free-form Fortran with FPM. It provides stationary and Krylov
iterative linear-system solvers for dense matrices and a package-local COO
sparse matrix type. No external BLAS or LAPACK installation is required.

## Implemented methods

Stationary methods:

- weighted Jacobi
- forward Gauss-Seidel
- SOR
- SSOR

Krylov methods:

- conjugate gradient (CG)
- biconjugate gradient (BiCG)
- BiCGSTAB
- conjugate gradient squared (CGS)
- Chebyshev semi-iteration
- restarted GMRES
- QMR

The `aux_fisch` generic generates the 2-D five-point discrete Poisson matrix in
either dense or COO sparse form.

## Build and test

```text
fpm build
fpm test
```

The source is standard free-form Fortran and has also been validated directly
with strict gfortran diagnostics and runtime checks. No system BLAS/LAPACK link
flags are used.

## Example

```fortran
use rlinsolve, only : dp, solver_result, lsolve_cg
real(dp) :: a(3,3), b(3), x0(3)
type(solver_result) :: out

! Fill A and b, then:
x0 = 0.0_dp
call lsolve_cg(a, b, out, xinit=x0)
```

`solver_result` contains `x`, `errors`, the upstream-style zero-based `iter`, a
status `flag`, and `converged`. Rank-2 right-hand-side matrices are also
supported; the corresponding `solver_matrix_result` stores all solutions plus
full per-column `solver_result` histories.

For several right-hand sides, use the same generic with a rank-2 `B`:

```fortran
use rlinsolve, only : dp, solver_matrix_result, lsolve_cg
real(dp) :: a(3,3), b(3,2)
type(solver_matrix_result) :: out

call lsolve_cg(a, b, out)
```

## Sparse representation

`sparse_matrix` stores one-based COO arrays `row`, `col`, and `val` together
with `nrow` and `ncol`. Duplicate COO entries are summed when converted to
dense form, and sparse matrix-vector products operate directly on the stored
entries.

## Translation coverage

Package status: **substantial**.

**12 of 12 (100%)** exported computational R functions are mapped. This
fraction measures mapped computational functions, not complete compatibility
with R classes, random-number behavior, warning/message text, sparse RHS Matrix
objects, or the RcppArmadillo sparse backend. Dense rank-2 right-hand-side
matrices are supported for every solver and are processed column by column, as
in the upstream wrappers.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `aux.fisch` | `rlinsolve` | `aux_fisch` | complete |
| `lsolve.bicg` | `rlinsolve` | `lsolve_bicg` | substantial |
| `lsolve.bicgstab` | `rlinsolve` | `lsolve_bicgstab` | substantial |
| `lsolve.cg` | `rlinsolve` | `lsolve_cg` | substantial |
| `lsolve.cgs` | `rlinsolve` | `lsolve_cgs` | substantial |
| `lsolve.cheby` | `rlinsolve` | `lsolve_cheby` | substantial |
| `lsolve.gmres` | `rlinsolve` | `lsolve_gmres` | substantial |
| `lsolve.gs` | `rlinsolve` | `lsolve_gs` | substantial |
| `lsolve.jacobi` | `rlinsolve` | `lsolve_jacobi` | substantial |
| `lsolve.qmr` | `rlinsolve` | `lsolve_qmr` | substantial |
| `lsolve.sor` | `rlinsolve` | `lsolve_sor` | substantial |
| `lsolve.ssor` | `rlinsolve` | `lsolve_ssor` | substantial |

See `API_COVERAGE.md` and `COMPATIBILITY.md` for details.

## Provenance and license

Rlinsolve is by Kisung You and is distributed under the MIT license. Original
package metadata and computational R/C++ sources are retained in `upstream/`
for provenance and are not compiled. See `NOTICE.md` and `LICENSE`.
