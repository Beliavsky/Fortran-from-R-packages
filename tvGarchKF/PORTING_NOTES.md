# Porting notes

## State-space recursion

The C routines `tvGarch1ab` and `tvGarch1ab2` differ only in whether the
intercept is treated as scalar or time-varying during the negative-state
repair. A one-coefficient polynomial produces a constant vector, so the two
paths are numerically equivalent in that case. The Fortran port uses the
array-valued intercept consistently.

The upstream sentinel value `100000000` for a missing squared return is
replaced internally by IEEE NaN detection. Appended prediction periods are
created as NaNs and receive the same no-measurement state update.

## Objective naming

The value called `loglik` upstream is the criterion

`sqrt(F_t) * r_t^2 / sigma_t^2 + log(sigma_t^2) - 0.5*log(F_t)`

summed over observed returns. It is minimized and is not the conventional
Gaussian log-likelihood. The Fortran result field is therefore named
`criterion`, while the compatibility routine retains `loglike` in its name.

## Constraints

Source-compatible mode checks only `alpha(t)+beta(t)<1`, as the R function
does in practice. `corrected_constraints=.true.` additionally requires
`omega(t)>0`, `alpha(t)>=0`, and `beta(t)>=0`. Corrected mode is recommended
for estimation and simulation.

## Deterministic functions

R accepts arbitrary expression strings through `parse` and `eval`. Embedding
an R expression interpreter is outside the numerical package. The two used
arguments (`u` and `3*(1-log(u))`) are built in, and arbitrary alternatives
can be supplied as pure scalar Fortran callbacks.

The upstream trigonometric routine applies the same trigonometric basis to
every coefficient after the intercept. This behavior is preserved; it does
not generate increasing harmonics automatically.

## Optimization

R chooses conjugate gradient for all-polynomial models and BFGS otherwise.
The supplied fGarch Fortran dependency provides a context-safe Nelder-Mead
optimizer, which is used here for all coefficient paths. Consequently,
parameter estimates need not be bit-for-bit identical to `stats::optim`.
The raw estimates and the upstream-style five-decimal rounded vector are both
returned.

## Rolling GARCH fits

`tv_parameter` calls the supplied `fit_garch11` implementation with zero mean
for the full sample and every overlapping block. Its output contains both the
local trajectories and the global stationary reference estimate; plotting is
left to the caller.
