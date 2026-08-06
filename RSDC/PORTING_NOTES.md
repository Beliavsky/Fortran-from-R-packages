# Porting notes

## Data layout

R matrices are translated to ordinary rank-2 Fortran arrays. Observations are
`y(T,K)`. Regime correlations are `rho(N,C)`, where `C=K(K-1)/2`. Correlation
vectors preserve R's `lower.tri` ordering, not row-major lower-triangle order.

## Transition timing

TVTP simulation and filtering both use the contemporaneous covariate row
`X(t,:)` to form `P_t`, matching RSDC 1.7-0. The first simulated state is drawn
from `transpose(P_1) * uniform_prior`. Smoothing uses `P_(t+1)`.

## Positive-definite search geometry

Natural pairwise correlations have a difficult joint constraint. Estimation
therefore searches canonical partial correlations in `(-0.999,0.999)` and maps
each point bijectively to a positive-definite correlation matrix. This follows
the same Joe/C-vine parameterization used in recent upstream global-search
code.

## Optimization

The supplied `deoptimr-modern-fortran` jDE solver is vendored under
`src/vendor`. It provides the stochastic global stage. A deterministic bounded
coordinate pattern search performs local refinement. This is an algorithmic
equivalent of the upstream DEoptim + L-BFGS-B workflow, not a bitwise port of
R's optimizers.

## Portfolio optimizers

For full-investment minimum variance, the unconstrained solution is
proportional to `inverse(Sigma) * 1`. Maximum-diversification weights are
proportional to `inverse(Sigma) * asset_volatility`. The long-only versions use
an active-set elimination method: negative-weight assets are removed and the
closed-form problem is resolved on the remaining set. This avoids external QP
and nonlinear-programming libraries.

## Inference

Scores and Hessians use central finite differences. OPG and sandwich matrices
use the same algebra as upstream. Parameter covariance can fail at boundaries
or singular optima; in that case `vcov` remains unallocated.

## License compatibility

RSDC declares GPL-3. The supplied DEoptimR translation is GPL-2.0-or-later and
is compatible with GPL-3.0-only distribution. The supplied mvtnorm translation
is GPL-2.0-only and is not linked or copied into the source graph. A small,
independently written Box-Muller/Cholesky multivariate-normal generator replaces
its only required role.
