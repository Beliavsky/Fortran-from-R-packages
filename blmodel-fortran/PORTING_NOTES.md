# Porting notes

## Data representation

The R input `dat = cbind(returns, probabilities)` is split into a rank-two
`returns` array and rank-one `probabilities` array. This avoids a special final
column and permits stronger dimension checking.

Result lists are represented by `moment_result`, `equilibrium_result`, and
`posterior_result` derived types.

## View-density callback

R's `match.fun` mechanism is replaced by a typed procedure callback. The
callback receives view points, the view mean, covariance, optional numeric
parameters, and returns allocated densities plus an integer status.

The R documentation says that points are stored in rows, but the implementation
uses `ncol(x)` as the number of points and subtracts a `k x n` matrix. The
Fortran API follows the implementation: points have shape `(k, n)` and each
column is one point.

## Risk inversions

The upstream implementation uses identical formulas for CVaR and deviation
CVaR, and identical formulas for LSAD and MAD. This behavior is retained.
The Fortran implementation handles empty tail sections explicitly rather than
using R's potentially fragile index ranges.

## Probability handling

Upstream routines assume scenario probabilities sum to one. The Fortran port
normalizes nonnegative probabilities with positive total mass. Therefore,
rescaling every prior probability by the same positive constant leaves results
unchanged.

## Linear algebra

Normal, Student-t, and power-exponential densities use Cholesky factorization.
This avoids explicit matrix inversion and determinant evaluation while
preserving the same formulas. Posterior view covariance matrices are explicitly
symmetrized before factorization.

## Frequency conventions

`returns_freq` is the number of observations per year. Sample covariance and
mean are annualized by multiplication. Equilibrium returns, covariance supplied
to the posterior, and user views are converted back to per-period values exactly
as in `BL_post_distr`.

## Licensing

The DESCRIPTION file declares `GNU General Public License version 3`, not a
later-version option. The translation therefore uses SPDX identifier
`GPL-3.0-only`.
