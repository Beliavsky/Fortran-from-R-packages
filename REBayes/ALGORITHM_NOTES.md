# Algorithm notes

## Finite-grid Kiefer-Wolfowitz solver

For the likelihood case used by almost all REBayes mixture routines, define

    B = A diag(d),  f >= 0,  sum(f) = 1.

The optimization problem is

    maximize  sum_i w_i log((B f)_i).

The port uses the exact EM multiplicative update

    f_j <- f_j * sum_i w_i B_ij / (B f)_i

with normalization and periodic vertex-direction line searches.  Convergence is
checked by the KW KKT condition

    max_j sum_i w_i B_ij / (B f)_i - 1 <= tol.

This removes the upstream dependency on MOSEK while solving the same finite-grid
NPMLE problem.  `kw_result%kkt_gap` exposes the stopping diagnostic.

## Noncentral Student t

For `Tncpmix`, the density is evaluated from the conditional-normal integral over
a chi-square variate.  A cached 96-point Gauss-Legendre rule on a transformed
semi-infinite interval is used.  This avoids a special-function library and is
deterministic.

## RLR

The regularized logistic problem

    min_theta L(theta) + lambda ||D theta||_1

is split with z = D theta and solved by ADMM.  The theta subproblem is solved by
damped Newton iteration; its Hessian is positive semidefinite plus the ADMM
quadratic term, and a small ridge is included for numerical rank deficiency.
The z step is exact soft thresholding.

## MEDDE

The upstream MOSEK formulation has the equality

    diag(dv) f + A z = e.

The Fortran port eliminates f:

    f = (e - A z) / dv,

and directly optimizes the same Shannon/log/power objective over the derivative
variables z, projected onto the upstream box constraints.  Backtracking rejects
steps that violate nonnegativity of f.  Dorder 1, 2, and 3 use exact finite-
difference coefficient matrices.

## Repeated-measures routines

R's `tapply` calls are replaced by `group_stats`, which computes group counts,
weighted means, weighted residual variances, total weights and sums of log
weights.  The likelihood formulas then follow the source routines directly.

## Interface differences

R histogram/collapse options are performance conveniences and are not required
for the likelihood.  The Fortran routines operate directly on observations and
optional weights.  Users can aggregate observations before calling the routines
if desired.
