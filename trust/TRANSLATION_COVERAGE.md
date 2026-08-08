# Translation coverage

## Translated directly

The computational content of `R/trust.R` is translated directly:

- trust-region iteration and radius update;
- minimization and maximization;
- optional parameter scaling;
- Hessian eigendecomposition;
- unconstrained Newton step when it lies inside the trust region;
- easy-easy, hard-easy, and hard-hard trust-region subproblem cases;
- one-dimensional secular-equation bracketing/root solve;
- actual/predicted change ratio `rho`;
- acceptance/rejection rules and radius quartering/doubling;
- `fterm` and `mterm` termination;
- restricted-domain trial points represented by signed infinity;
- iteration/path diagnostics corresponding to `blather=TRUE`.

## Fortran-specific replacements

R's `eigen(..., symmetric=TRUE)` is replaced by a self-contained cyclic Jacobi symmetric eigensolver.  This avoids a mandatory LAPACK dependency while retaining the eigendecomposition algorithmic structure used by upstream `trust`.

R's `uniroot` is replaced by a safeguarded bisection solve of the same monotone secular equation.  Upstream already brackets the root analytically; bisection therefore preserves the trust-subproblem solution without depending on R's root-finder implementation.

R exceptions/`try()` are represented by the objective callback's integer `status`.  A nonzero status causes a controlled unsuccessful return rather than a language-level exception.

R lists and `blather` matrices/vectors are represented by Fortran derived types with allocatable components.

## Omitted

There is no plotting code in the upstream package.  R-specific list checking, condition/warning machinery, and `...` argument forwarding are not reproduced; typed Fortran procedure interfaces replace them.

## Numerical equivalence

The algorithm and decision thresholds follow the upstream R source.  Floating-point trajectories can differ slightly because the Fortran release uses its own symmetric eigensolver and bisection implementation rather than R's LAPACK-backed `eigen` and `uniroot`.
