# Porting notes

## Fortran representation

R vectors and matrices are represented by rank-one and rank-two `real(dp)`
arrays. R lists returned by fitting procedures are replaced by derived types
with explicit status fields. R functions passed as quantile or CDF arguments
are represented by typed procedure callbacks.

## Numerical libraries

The port is dependency-free. It includes native implementations of the
probability functions, matrix solves, Cholesky decomposition, numerical
Hessians, bisection, and Nelder-Mead optimization needed by the package.

## Corrections and safeguards

- The tied-data branch of `mean_excess_np` is written without assuming
  short-circuit evaluation. This prevents an invalid `z(0)` reference under
  checked Fortran execution.
- Composite-distribution interval boundaries are evaluated through explicit
  conditionals rather than `merge`, because Fortran may evaluate both
  arguments.
- Block rearrangement uses a temporary block matrix, avoiding in-place indexed
  assignments whose result can depend on evaluation order.
- Rearrangement routines avoid R's `Inf` endpoint quantiles by using the same
  midpoint replacements documented upstream.
- Covariance and Hessian inverses are explicitly symmetrized.
- Brownian routines validate time grids and dimensions and accept deterministic
  uniforms or seeds.
- Black-Scholes handles zero time and zero volatility directly.
- GARCH fitting validates positive conditional variance and standardized-t
  degrees of freedom.

## Formula compatibility

- The GARCH likelihood follows the upstream Zumbach reparameterization and
  method-of-moments fixed unconditional variance.
- The default GARCH normal likelihood includes the full normalizing constant,
  as does the R implementation.
- Sample quantiles use the discrete type-1 convention used by the translated
  nonparametric risk routines.
- `mean_excess_np` intentionally follows the upstream unique-value calculation
  and then expands tied thresholds according to their observed counts.
- The adaptive rearrangement interface uses a margin-indexed quantile callback
  instead of an R list of functions.

## Concurrency

Objective data for native fitters are stored in module scope to avoid nested
procedure trampolines and executable-stack requirements. Calls to the same
fitter should therefore be serialized when used from multiple threads.
