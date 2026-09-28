# Compatibility audit

This file distinguishes three targets that can otherwise be conflated when
translating sanic 0.0.2: the documented R API, the algorithms intended by its
C++ branches, and the literal behavior of the released source. The Fortran
translation targets the documented computational behavior unless a literal
upstream quirk is both intentional and useful.

## Upstream defects intentionally corrected

Several C++ files declare a solver object before an `if` chain, then redeclare a
different solver with the same name inside a branch. Standard C++ block scope
means the inner declaration is destroyed at the end of the branch and the
outer default solver is the object that is subsequently `compute()`d and
`solve()`d. Consequently, the released code does not activate several choices
that the R documentation explicitly advertises.

| Upstream source | Documented choice | Literal sanic 0.0.2 effect | Fortran behavior |
|---|---|---|---|
| `upstream/src/solve_Cholesky.cpp` | dense LLT vs LDLT; sparse LLT/LDLT and ordering | branch-local objects shadow the outer default solver | implements LLT/LDLT and sparse ordering choices |
| `upstream/src/solve_LU.cpp` | full-pivot dense LU; sparse AMD/COLAMD/natural ordering | branch-local solver is discarded | implements full pivoting and each ordering selector |
| `upstream/src/solve_QR.cpp` | unpivoted/column/full-pivot dense QR; sparse ordering | branch-local solver is discarded | implements distinct dense selectors and sparse ordering |
| `upstream/src/solve_iterative.cpp` | identity/diagonal/incomplete preconditioners | branch-local iterative solver is discarded | implements the documented preconditioner choices |
| `upstream/src/svd.cpp` | alternate Jacobi QR preconditioners | branch-local SVD object is discarded | public R API does not expose this selector; Fortran uses the documented default path |

The R wrapper `upstream/R/solve_iterative.R` contains two additional issues.
When `iter` is supplied it calls `int_check(tol, ...)` instead of
`int_check(iter, ...)`, so the requested iteration count is not honored. Its
`dim(a)[1] != length(b)` check also rejects a matrix RHS with more than one
column even though the underlying Eigen routines accept multiple RHS columns.
Fortran honors `iter` and supports multiple RHS columns.

The generic `solve2()` path has a related sparse-wrapper defect: it calls
`maybe_symmetric()` before dispatch, and that helper explicitly rejects inputs
that are not base R matrices. Thus a sparse `dgCMatrix` can fail before reaching
the sparse direct solvers even though the underlying solver functions support
sparse input. `solve2_sparse` implements the documented sparse dispatch rather
than preserving that wrapper failure.

These are deliberate bug fixes, not accidental incompatibilities. Reproducing
them would make the translated API contradict sanic's own documentation.

## Eigen semantics matched explicitly

- Iterative solvers use Eigen's default maximum-iteration convention of twice
  the number of matrix columns when the caller supplies no positive `iter`.
- Explicit self-adjoint eigensolves use only the input diagonal and lower
  triangle, reflecting the lower triangle to form the symmetric matrix. This
  matches `Eigen::SelfAdjointEigenSolver`, which documents that only the lower
  triangular part is referenced.
- Rank-deficient/underdetermined pivoted QR returns a valid basic/arbitrary
  solution rather than silently substituting a pseudoinverse solution.

## Remaining numerical differences

- Dynamic minimum-degree orderings are package-local approximations, not exact
  Eigen/SuiteSparse AMD or COLAMD permutations.
- Package-local ILUT and incomplete Cholesky follow Eigen's main drop/fill,
  scaling, permutation, and shifting policies but are not guaranteed to make
  identical numerical choices.
- The BDC SVD selector uses bidiagonalization plus implicit QR, not Eigen's
  divide-and-conquer `BDCSVD`.
- The general eigensolver is a portable Hessenberg/shifted-QR implementation,
  so eigenvalue ordering and roundoff can differ from Eigen.
- Missing Arnoldi/Lanczos start vectors are deterministic in Fortran instead of
  consuming R's random-number stream. Early-breakdown inactive outputs are
  zero-filled rather than leaving Eigen dynamic-vector tail entries undefined.
- R list/class construction, warnings, argument partial matching, and printed
  diagnostics are not reproduced. Numerical status is returned through Fortran
  arguments and result types.

## Validation policy

Regression tests emphasize invariants that are meaningful across Eigen and
Fortran implementations: solve residuals, reconstruction residuals, sorted
spectra where ordering is not contractual, requested output shapes, explicit
selector behavior, and convergence diagnostics. Cases that encode a known
upstream defect are documented as deliberate divergences instead of being
used as golden outputs.
