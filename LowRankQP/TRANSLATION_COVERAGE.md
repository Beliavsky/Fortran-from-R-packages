# Translation coverage

## Fully translated computational surface

LowRankQP 1.0.6 has one public numerical routine, `LowRankQP()`. Its complete
numerical C core is represented in the Fortran package:

- initialization of strictly interior primal/dual variables;
- primal and dual residuals;
- complementarity, duality gap, and termination statistic;
- predictor and corrector directions;
- equality-constraint Schur complement solve;
- fraction-to-boundary step rule;
- LU path;
- Cholesky path;
- SMW path;
- product-form Cholesky factorization and solve;
- upper/lower box constraints;
- optional absence of equality constraints.

R dimension checks are represented by Fortran result status codes.

## Deliberately omitted

- `.C` registration and R native-routine registration;
- R list construction;
- `Rprintf` formatting code as an API requirement (Fortran `verbose` output is available);
- BLAS/LAPACK wrappers whose only purpose was to expose operations already
  represented directly in modern Fortran.

There is no plotting code in the package.

## Source behavior retained

The upstream source uses a special convention:

- square `V`: the Hessian is `H = V`;
- non-square `V`: the Hessian is `H = V V^T`.

This is retained exactly. Consequently the methods recommended by the manual
are `LU`/`CHOL` for square `V` and `SMW`/`PFCF` for low-rank non-square `V`.
The R default is nevertheless `PFCF`; the Fortran default retains that choice
for source compatibility rather than silently changing it.

The original C source uses BLAS/LAPACK LU and Cholesky. The Fortran port uses
self-contained dense factorizations with the same mathematical operations;
roundoff trajectories can therefore differ slightly while converging to the
same QP solution.
