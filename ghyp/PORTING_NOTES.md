# Porting notes

## Parameterization

The internal representation follows the normal mean-variance mixture

```text
W ~ GIG(lambda, chi, psi)
X = mu + W gamma + sqrt(W) A Z
A A' = scatter
Z ~ N(0, I)
```

For a univariate chi/psi constructor, the user-facing `sigma` argument is the
standard-deviation scale, so the stored scatter value is `sigma**2`.

## Density formulas

The implementation uses closed-form generalized-hyperbolic densities based on
real-order modified Bessel K functions. Student and variance-gamma boundary
cases are evaluated with their limiting formulas. Symmetric Student t uses a
closed-form t density and CDF.

## Native Bessel K

Fortran does not provide a standard real-order modified Bessel K intrinsic.
The port evaluates

```text
K_nu(x) = integral_0^infinity exp(-x cosh(t)) cosh(nu t) dt
```

using cached Gauss-Legendre nodes and a scaled integrand. The half-order closed
form is used when applicable.

## GIG simulation

The upstream package contains a specialized C rejection sampler. The port uses:

- direct gamma generation when `chi = 0`;
- direct inverse-gamma generation when `psi = 0`;
- inverse-CDF generation for a general GIG law.

The general method is slower but has a compact, dependency-free implementation.

## Multivariate probabilities

The upstream `pghyp` function estimates multivariate probabilities by drawing
from the fitted law. `pghyp_rectangle` preserves that methodology and also
returns its binomial Monte Carlo standard error.

## Fitting

The upstream multivariate fitter uses family-specific EM calculations and R
optimization helpers. This port uses one typed likelihood representation and a
native Nelder-Mead optimizer for univariate and multivariate models.

Positive parameters and Cholesky diagonal entries are optimized on logarithmic
scales. Student degrees of freedom are represented as

```text
nu = 2 + exp(theta)
```

which ensures finite variance during fitting.

The optimizer objective data are held in module state. Concurrent fitting calls
must be serialized.

## Risk integration

Expected shortfall, Omega, and noninteger or absolute moments use deterministic
quadrature over transformed infinite intervals. Integer moments are evaluated
analytically from GIG raw moments and Gaussian conditional moments.

## Portfolio optimization

When risk is standard deviation, or when the distribution is symmetric, the
upstream analytical mean-variance structure is preserved. Non-symmetric VaR
and ES objectives use unconstrained weight parameterization with the first
weight enforcing the full-investment condition.

Target-return VaR/ES optimization uses a quadratic return penalty rather than
the upstream algebraic elimination of two weights. SD target-return problems
continue to use an exact KKT system.

## Thread safety

Distribution evaluation and simulation are reentrant after quadrature-rule
initialization. The fitting and non-symmetric portfolio optimizers use
module-held objective context and should be externally serialized.

## Numerical safeguards

- Scatter matrices are Cholesky-validated.
- Covariance and inverse-Hessian outputs are explicitly symmetrized.
- Quantile brackets expand automatically.
- Tail integrations avoid evaluating infinite endpoints.
- Variance-gamma singularities at the location are handled by their analytic
  limit when finite.
- Invalid parameter combinations return typed error states rather than relying
  on R exceptions.
