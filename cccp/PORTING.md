# Porting notes

## Scope

The numerical behavior represented by all exported R functions is included:

| R function or method | Fortran API |
|---|---|
| `ctrl` | `ctrl`, `cccp_control` |
| `nnoc` | `nnoc` |
| `socc` | `socc` |
| `psdc` | `psdc` |
| `nlfc` | `nlfc` |
| `dlp` | `dlp`, `solve_lp` |
| `dqp` | `dqp`, `solve_qp` |
| `dnl` | `dnl`, `solve_dnl` |
| `dcp` | `dcp`, `solve_dcp` |
| `cccp` | generic `cccp` from module `cccp_api`, or `cccp_solve` from module `cccp`; `dnl`/`dcp` for callbacks |
| `gp` | `gp` |
| `l1` | `l1` |
| `rp` | `rp` |
| `cps` | generic `cps` for problem objects |
| extractor generics | `getx`, `gety`, `gets`, `getz`, `getstate`, `getstatus`, `getniter`, `getparams` |

R reference classes, S4 dispatch, Rcpp modules, formatted `show` methods, and R
list/matrix coercion are presentation/runtime infrastructure and are not ported.

## Solver substitution

The original C++ implementation is partially derived from CVXOPT and uses a
primal-dual predictor-corrector interior-point method with Nesterov-Todd scaling
and a homogeneous-style set of primal/dual diagnostics.

This translation uses a different but mathematically equivalent problem-level
approach:

1. phase-I strict-interior recovery;
2. primal logarithmic barriers for every cone and nonlinear inequality;
3. equality-constrained Newton steps through a dense KKT system;
4. backtracking line search; and
5. increasing barrier weights until the theoretical barrier gap is below the
   requested tolerance.

This substitution avoids Rcpp, RcppArmadillo, CVXOPT-derived object machinery,
and external optimization packages. BLAS/LAPACK are used only for dense linear
solves, Cholesky inversion, and symmetric eigenvalues.

Consequences:

- primal solutions should agree to numerical tolerance;
- iteration counts do not match the R package;
- `y` and `z` are barrier-based approximate multipliers;
- infeasibility certificates are less comprehensive than the original
  primal-dual method;
- problems without a strict interior can be harder than under a homogeneous
  self-dual embedding; and
- dense storage is intended for small and medium problems, matching the package's
  examples rather than large sparse conic programs.

## Cone representation

The sign conventions are unchanged:

- NNOC: `h - G*x` must be elementwise nonnegative.
- SOCC: the first component of `h - G*x` is the Lorentz scalar.
- PSDC: `h - G*x` is the column-major vectorization of a symmetric PSD matrix.

Fortran arrays are column-major, matching R's default matrix storage and the
original Armadillo implementation.

## Nonlinear callbacks

R lists of scalar function, gradient, and Hessian closures are replaced by one
batch callback that fills all constraint values, gradients, and Hessians. This
removes runtime type dispatch and substantially reduces callback overhead.

An initial point must be in the callback's mathematical domain. It need not be
strictly feasible: the port augments nonlinear constraints with the same phase-I
slack used for cone constraints.

## Geometric programming

The stable log-sum-exp formulas and analytical derivatives from the original
`fgp` routine are preserved. The Fortran wrapper uses temporary module context
for callback dispatch; therefore simultaneous threaded calls to `gp` or `rp`
from the same process are not reentrant. Independent processes are unaffected.

## Risk parity

The original wrapper doubles `P` before evaluating

```text
0.5*x^T*P*x - sum(mrc*log(x)).
```

The Fortran implementation preserves that convention and then normalizes the
positive optimizer to unit budget.
