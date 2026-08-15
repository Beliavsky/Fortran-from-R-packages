# Translation notes

`LaplacesDemon` is callback-centric: the R package expects `Model(theta, Data)`
to return an R list containing the log posterior, monitored quantities,
deviance and other bookkeeping. The Fortran port exposes the numerical
contract directly through explicit procedures and arrays, for example:

```fortran
function log_target(theta) result(lp)
  real(dp), intent(in) :: theta(:)
  real(dp) :: lp
end function
```

Data may be held in a module or otherwise made available to the callback. This
avoids recreating dynamic R lists while preserving the numerical algorithms.

## Derivatives

Finite-difference gradients, Hessians and Jacobians are used where the
translated algorithm needs derivatives and no analytic callback is supplied.
This includes MALA, HMC/HMCDA, NUTS, several Laplace-approximation optimizers
and Salimans2 variational Bayes.

BHHH differs from the scalar-target optimizers because its update is defined
from per-observation likelihood components. The Fortran routine therefore
takes an explicit callback that returns those components rather than trying to
recover them from an R model object.

## MCMC adaptation and callback differences

NUTS follows the upstream Hoffman-Gelman structure: a slice variable in joint
position/momentum density, recursive binary-tree doubling, reservoir selection
among valid states, U-turn stopping and dual averaging. HMCDA uses the upstream
dual-averaging constants.

DRAM preserves the upstream two-stage proposal and adaptive covariance update.
RAM follows the Vihola-style covariance-factor adaptation. `pcn_sample`
deliberately follows LaplacesDemon's callback-difference acceptance convention;
callers wanting prior-reversible pCN should provide the corresponding
likelihood component as the target callback.

The t-walk retains the upstream two-state construction and traverse, walk,
blow and hop kernels. RDMH uses the upstream random-dive proposal
`epsilon = U**s`, with `U` uniform on `(-1,1)`, `s` in `{-1,+1}`, and the
associated Jacobian.

Gibbs sampling takes an explicit conditional-draw callback. INCA is represented
by an in-process multiple-chain implementation; the R package's external
network/HPC exchange mechanism is orchestration rather than a separate
probability kernel. SAMWG/SMWG updating variants similarly take the baseline
and dynamic-index information explicitly.

The RJ routine implements the package's fixed-vector variable-selection
semantics with an active mask, parameter-inclusion probabilities and a model
size prior. It does not attempt to emulate R's ability to mutate arbitrary
object structures or Fortran array rank dynamically.

## Quadrature

Three complementary approximation paths are available:

- tensor iterative Gauss-Hermite, retained with a conservative dimension cap
  because tensor node counts grow exponentially;
- componentwise iterative quadrature;
- adaptive sparse-grid quadrature, currently capped at eight dimensions.

The latter two were added in v0.3 and replace the earlier v0.2 note that they
were future work.

## Distribution compatibility

The distribution API is scalar/matrix-first rather than R-vector-recycling
first. Where upstream offers covariance, precision and Cholesky forms, the
Fortran implementation uses explicit SPD matrix transformations. Exported
`distributions.R` kernels are represented; local nested helpers are folded
into their caller rather than exposed solely to mirror an internal R symbol.

A few formulas in the upstream R helpers are numerically fragile. The Fortran
port uses stable log-sum-exp/log1p-style evaluation where this does not change
the mathematical distribution. Generic truncation moments are computed by
numerical quadrature over the supplied scalar PDF/CDF callbacks.

## Deliberately omitted infrastructure

Plotting, printing, R classes/lists, formula/model-frame handling, monitor and
file/database storage, package loading, environment mutation and cluster/socket
transport are intentionally not reproduced. These are application/runtime
services rather than missing numerical algorithms.
