# Porting notes

## Naming and result representation

R dots were converted to underscores only where needed. The exported operator
names `L`, `A`, `D`, `Lstar`, `Astar`, and `Dstar` are valid Fortran identifiers
and were retained. R lists are represented by `type(graph_result)`.

## Replaced dependencies

- Armadillo/Eigen symmetric eigensolvers were replaced by a Jacobi method.
- `MASS::ginv` was replaced by a symmetric eigenvalue pseudoinverse.
- `quadprog::solve.QP` was replaced by either an exact Euclidean simplex
  projection or a projected-gradient nonnegative quadratic solver, according to
  the original problem structure.
- `stats::isoreg` was replaced by pool-adjacent-violators isotonic regression.
- `CVXR` in `learn_graph_sigrep` was replaced by a Laplacian edge-weight
  parameterization and projected convex quadratic optimization.
- `qr.Q(qr(1), complete=TRUE)[,2:p]` was replaced by a Helmert basis for the
  orthogonal complement of the all-ones vector.

## Deliberate numerical details

The signal-representation Laplacian update enforces the original trace and
Laplacian constraints through nonnegative edge weights with a negligible
positive edge floor. This is slightly stronger than requiring only positive
node degrees, but avoids an external general-purpose conic solver.

`learn_smooth_graph` clips its final extragradient edge update at zero. The
underlying method is constrained to nonnegative graph weights; the clipping
prevents tiny finite-precision violations.

The original R `accuracy()` wrapper returns the specificity element by mistake.
The Fortran `accuracy()` returns the actual accuracy element computed by
`metrics`.

## Dense implementation

The original package uses sparse-capable R libraries and C++ linear algebra.
This translation is intentionally dependency-free and dense. It is appropriate
for validation, research prototypes, and moderate graph sizes; large sparse
problems should eventually use optimized BLAS/LAPACK and sparse matrix backends.
