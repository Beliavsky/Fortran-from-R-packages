# Translation coverage

## Directly covered

The complete exported computational surface of limSolve is represented at
array level: `Solve`, `Solve.banded`, `Solve.block`, `Solve.tridiag`, `ldei`,
`ldp`, `linp`, `lsei`, `nnls`, `resolution`, `varranges`, `varsample`,
`xranges`, and `xsample`.

## Architectural differences

1. `nnls` is a modern Lawson-Hanson-style active-set implementation.  It does
   not retain the original in-place Householder storage layout of `xNNLS`.
2. `ldp` and the reduced constrained step in `lsei` are expressed as convex
   QPs and solved with the vendored Goldfarb-Idnani quadprog translation.
   The R package's type-1 path instead calls the older Lawson-Hanson `xLDP` /
   `xDLSEI` code directly.  Equality-inconsistent LSEI problems still use a
   generalized-inverse equality solution before optimization.
3. `tolrank` is represented by the common numerical rank tolerance rather than
   two separately configurable Lawson-Hanson rank thresholds.
4. `fulloutput` covariance is computed from the reduced normal matrix rather
   than reproducing `xDLSEI`'s packed work-array covariance output bit for bit.
5. `Solve` and `resolution` use a self-contained symmetric-Jacobi
   pseudoinverse/SVD-equivalent construction instead of MASS `ginv`/R `svd`.
6. `Solve.banded` and `Solve.tridiag` reconstruct/solve with modern pivoted
   dense elimination.  `Solve.block` reconstructs the almost-block-diagonal
   matrix and solves it densely rather than preserving COLROW's no-fill
   alternating row/column factorization.  Results are equivalent but memory
   scaling differs for large systems.
7. `xsample` implements all three walk kernels and the null-space likelihood
   construction.  The first release does not reproduce the R wrapper's hidden
   equality discovery pass or its optional posterior gamma sampling of an
   unknown data standard deviation.  Automatic jump selection is simplified.
8. R data frames, column names, warnings, printing, and plotting/vignette code
   are intentionally omitted.

The full original source, including Lawson-Hanson, LINPACK/LAPACK, and COLROW
routines, is retained in `original/limSolve-master` for provenance and future
line-for-line modernization work.
