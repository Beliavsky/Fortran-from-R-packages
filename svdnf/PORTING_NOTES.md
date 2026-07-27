# Porting Notes

## Array orientation

The R implementation stores volatility nodes in ascending order but reverses
filter vectors in several expressions. The Fortran implementation stores and
propagates probabilities in natural ascending-node order. This removes repeated
`N:1` indexing and makes transition matrices column-stochastic.

## Transition probabilities

The DNF recursion uses normal interval probabilities for the latent-state grid,
matching `probCalculator`. Prediction percentiles also use these interval
probabilities. The upstream percentile helper instead constructs a matrix from
point densities and normalizes afterward; using the same discretized transition
law as the filter is more internally consistent.

## Forecast starting distribution

Upstream `predict.SVDNF` samples from filter column `T`, even though the final
posterior is stored in column `T + 1`. `predict_filter` uses the final posterior
by default. Set `upstream_filter_index=.true.` to reproduce the earlier column.

## Custom models

R discovers custom parameters by inspecting function argument names. Fortran
cannot safely reproduce run-time argument introspection, so custom dynamics use
explicit parameter arrays and typed procedure callbacks. Custom estimation uses
a caller-supplied parameter-setter callback.

## Optimization

R delegates to `stats::optim`, whose method is selected by caller arguments. The
Fortran port supplies a self-contained Nelder-Mead optimizer and numerical
Hessian. Random streams and convergence paths therefore need not match R
exactly, while likelihood evaluation is based on the translated DNF equations.

## Stationarity constraints

The R constructor allows `phi` in `(-1,1)`, but the R optimizer rejects
non-positive `phi`. The Fortran validation follows the stationary interval
`abs(phi) < 1`, allowing negatively autocorrelated log volatility.

## Jump-size discretization

For zero jumps, the gamma jump-size distribution is treated as a point mass at
zero. This avoids passing gamma shape zero to a CDF, which can produce undefined
values in generic probability libraries.

## Underflow handling

If a likelihood contribution underflows, the result is marked invalid and its
log likelihood is set to negative huge. This is preferable to returning a
partially updated filter as a nominally successful result.

## Simulation

The model equations and jump conventions are preserved. Optional seeds provide
repeatability within the Fortran implementation; they do not reproduce R's RNG
streams.

## Factor matrices

Fortran factor matrices are shaped `(observations, factors)`. Forecast factor
matrices are shaped `(horizon, factors)`. This is explicit and avoids the
orientation guessing used by R methods.

## Thread safety

Built-in and custom optimization use module-held objective context to avoid
nested-procedure trampolines and executable-stack linker warnings. Concurrent
optimizer calls from multiple threads should therefore be externally
serialized. Filtering and simulation calls themselves do not use this context.
