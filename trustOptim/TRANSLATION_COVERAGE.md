# Translation coverage

## Directly translated computational behavior

- Steihaug (1983) truncated/preconditioned conjugate-gradient trust-region
  subproblem solver.
- Negative-curvature termination at the trust-region boundary.
- Trust-boundary intersection calculation in the preconditioner-induced norm.
- Actual/predicted improvement ratio.
- Contract/move/expand acceptance states and radius updates.
- Gradient-norm, minimum-radius, and maximum-iteration stopping rules.
- SR1 Hessian updates and the upstream `1e-7` denominator acceptance test.
- BFGS Hessian updates.
- Identity and full-Cholesky quasi-Newton preconditioning.
- User-supplied sparse Hessian mode.
- Upstream scaled modified-Cholesky construction for sparse preconditioning.
- `function.scale.factor` minimization/maximization semantics.
- The exported binary-choice value, gradient, and sparse Hessian routines,
  including both `order.row` layouts.

## Architectural changes

### R/Rcpp/Eigen interface

R `Function`, `NumericVector`, `dgCMatrix`, Rcpp registration, interrupt
checking, and S4 conversion are replaced by explicit Fortran procedure
interfaces, allocatable arrays, and a native sparse symmetric type.

### Sparse factorization

Upstream uses Eigen `SimplicialLLT` for the sparse modified-Cholesky
preconditioner.  This translation keeps the Hessian sparse during all CG
Hessian-vector products, but converts the scaled preconditioner matrix to dense
form for a self-contained Cholesky factorization.  Therefore the optimizer is
algorithmically sparse in its trust-region operator but the optional
preconditioner can require O(n^2) memory.

### Cholesky failure handling

Eigen's quasi-BFGS Cholesky path does not explicitly check factorization
failure in the upstream code.  The Fortran port falls back to identity
preconditioning if the dense quasi-Hessian is not positive definite.  This is
a defensive portability behavior, not a change to the trust-region model.

### BFGS division guards

The upstream BFGS formulas divide directly by `y^T s` and `s^T B s`.  The
Fortran port skips the update if either denominator is at machine-underflow
scale, avoiding NaNs without changing ordinary trajectories.

### Reporting

Formatted R console reports (`report.level`, column widths, periodic headers)
are not reproduced.  The controls remain in the type for API correspondence,
and the final result includes the last CG iteration count and stop reason.
This is presentation/I/O rather than numerical behavior.

## Omitted R infrastructure

- Rcpp/RcppEigen registration and object conversion.
- Matrix S4 class conversion.
- R user-interrupt handling.
- Package documentation website/vignettes as executable R documents.

The complete original source tree is retained under `original/`.
