# NOTICE and provenance

This directory is a modern Fortran translation of the computational code in the
R package **sanic**, version 0.0.2 (published on CRAN 2023-08-23).

Upstream project: <https://github.com/nk027/sanic>

Upstream authors and contributors recorded in `DESCRIPTION`:

- Nikolas Kuschnig, author and maintainer
- Lukas Vashold, contributor
- Yixuan Qiu, contributor

The upstream package is distributed under GPL (>= 3), with copyright (C) 2023
Nikolas Kuschnig. The original `LICENSE` file is retained verbatim at the root
of this translation. A copy of GPLv3 is included as `LICENSE.GPL-3` when
available from the build environment.

The upstream package also carries and uses `UpperHessenbergEigen.h` from the
Spectra C++ library, with copyright notices for Gael Guennebaud, Jitse Niesen,
and Yixuan Qiu under MPL-2.0. The exact upstream header is preserved at
`upstream/src/UpperHessenbergEigen.h` for provenance only and is **not compiled
or used** by the Fortran implementation. The Arnoldi projected-eigenvalue
operation is independently implemented in Fortran. The original sanic notice
is retained in `LICENSE`, and `LICENSE.MPL-2.0` contains the MPL-2.0 text.

For provenance and coverage review, `upstream/` contains the upstream R files,
the sanic-owned C++ translation targets, and the preserved Spectra header used
by the original package. Generated Rcpp glue remains excluded as R-specific
interface code.

## Dependency review

Before implementation, the `Fortran-from-R-packages` repository was checked for
shared linear-algebra facilities. The repository contains `rfortran-linalg`,
which is appropriate for many general dense linear-algebra calls. This
translation does not add a dependency in its present revision because sanic's
public computational surface is itself a collection of specifically selectable
LLT/LDLT/LU/QR/Krylov algorithms and sparse matrix operations. The package-local
implementations keep those choices visible and make the current translation
self-contained. A future backend revision can replace common dense kernels with
`rfortran-linalg` where the shared API reproduces the requested sanic semantics
without hiding pivot/method differences.

No BLAS, LAPACK, ARPACK, Rcpp, RcppEigen, Eigen, Spectra, `r.f90`, `r_mod.f90`,
or `rfortran-compat` source is copied into the translated package. The FPM
manifest has no system-library links.

## Material compatibility differences

- Sparse Cholesky and LU factor dynamic sparse rows and permit fill-in. Their
  dynamic minimum-degree graph orderings approximate, but do not reproduce,
  Eigen/SuiteSparse AMD and COLAMD exactly. Natural ordering is supported
  directly, and the ordering graph itself uses sparse adjacency lists.
- Sparse QR remains sparse during the numerical factorization, using Givens
  rotations, structural ordering, and numerical column pivoting. Wide systems
  return a valid rank-revealing basic solution rather than forcing a
  pseudoinverse/minimum-norm solution. Its structural ordering still differs
  from Eigen's exact AMD/COLAMD implementations.
- `solve_cg(..., precond=2)` uses package-local dual-threshold ILUT with fill
  limiting for BiCGSTAB and scaled, permuted, adaptively shifted limited-memory
  incomplete Cholesky for CG. These follow the main Eigen algorithmic policies
  but are not claimed numerically identical to Eigen's `IncompleteLUT` and
  `IncompleteCholesky`. LSCG follows upstream behavior and keeps its default
  diagonal preconditioner when `precond=2` is requested.
- General eigendecomposition reduces to upper Hessenberg form and uses shifted
  complex Householder QR plus inverse iteration rather than Eigen's
  `EigenSolver`. Explicit self-adjoint solves mirror Eigen's documented rule
  that only the lower triangular part of the input defines the matrix.
- The two SVD selectors use distinct portable paths. `BDC` uses Householder
  bidiagonalization followed by implicit bidiagonal QR iteration. `Jacobi` uses
  column-pivoted Householder QR preconditioning and two-sided Jacobi sweeps,
  with a one-sided Jacobi convergence fallback. The `BDC` path is not Eigen's
  divide-and-conquer implementation.
- When Arnoldi/Lanczos omits the starting vector, Fortran uses a deterministic
  vector proportional to `1:n` rather than R's `rnorm()` stream.
- The sanic 0.0.2 C++ sources contain block-scope solver declarations that
  shadow the solver object used after the branch. As a result, several literal
  upstream pivot, ordering, and preconditioner selectors are accidentally
  ignored. This translation follows the documented intended selector behavior
  rather than reproducing those shadowing defects. The R `solve_cg()` wrapper
  also validates an explicit `iter` from `tol` by mistake; Fortran honors the
  caller's `iter` value. See `COMPATIBILITY.md`.
- Plotting, R classes, R warnings, verbose console presentation, Rcpp glue, and
  other R-specific interface behavior are outside the computational translation.
