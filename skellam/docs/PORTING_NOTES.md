# Porting notes

## Numerical implementation

The exact PMF evaluates the integer-order modified Bessel function of the first
kind with a positive power series accumulated in log space. A fixed-order
asymptotic expansion is reserved for very large arguments.

The R package computes the CDF through noncentral chi-square identities and
falls back to a saddlepoint approximation outside the working range. The
Fortran exact CDF instead builds a conservative finite support, evaluates the
log-PMF across that support, and normalizes with log-sum-exp. This removes the
need for a noncentral chi-square implementation and gives stable lower, upper,
and log tails over the tested range.

Quantiles are obtained from the same normalized support. The support width is
expanded according to the requested tail probability.

Poisson random generation uses inversion for small rates and the PTRS
transformed-rejection algorithm for larger rates.

## Estimation

The two-rate MLE preserves the upstream reduction

```text
lambda1 = lambda2 + sample_mean
```

and optimizes a log-transformed free parameter. BFGS is used first with a
Nelder-Mead fallback. Standard errors are calculated from the unconstrained
log-rate Hessian at the optimum.

Regression uses separate exponential links for the two latent Poisson rates.
The default response ordering follows the upstream regression example:
`response = count2 - count1`. A direct `count1 - count2` convention is available
through an optional argument.

## R behavior not reproduced

- Vector recycling of unequal-length arguments
- R warnings for noninteger PMF inputs
- Formula and data-frame conversion
- `options(verbose=TRUE)` validation behavior
- R-specific `NaN`, `Inf`, names, and class/list presentation

The PMF API therefore accepts integer quantiles explicitly. Invalid inputs are
reported through optional integer status arguments where available.
