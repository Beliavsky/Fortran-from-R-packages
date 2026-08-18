# Porting notes

## Distribution

A generalized Hermite variable has the exact representation

`X = P1 + m*P2`

with independent `P1 ~ Poisson(a)` and `P2 ~ Poisson(b)`. Consequently

- mean = `a + m*b`
- variance = `a + m^2*b`.

The upstream package evaluates small-parameter probabilities through a
recurrence and switches to Edgeworth/Cornish-Fisher approximations when either
`a` or `b` exceeds 20. Both paths are translated:

- `int_hermite` preserves the recurrence;
- `edg` preserves the Edgeworth CDF approximation;
- `cofi` preserves the Cornish-Fisher approximate quantile.

For validation and robust fallback, the Fortran library also exposes an exact
log-sum PMF (`dhermite_exact`) based on the independent-Poisson representation.

`dhermite`, `phermite`, and `qhermite` preserve the upstream threshold by
default. Pass `exact=.true.` to force exact probabilities/quantiles.

### Upstream fallback correction

When the upstream Edgeworth PMF difference becomes nonpositive, `dhermite`
falls back to `pnorm(x+0.5)-pnorm(x-0.5)` with standard-normal coordinates,
ignoring the Hermite mean and variance. The Fortran translation instead falls
back to the exact Hermite PMF. This affects only that pathological fallback.

The public Edgeworth CDF result is clipped to `[0,1]`. Quantile endpoints are
handled explicitly (`0 -> 0`, `1 -> huge(int64)`).

## Regression parameterization

`glm.hermite` parameterizes the distribution by fitted mean `mu`, integer order
`m`, and dispersion index `d = Var(Y)/E(Y)`.

For `m > 1`, the equivalent Poisson-component parameters are

`b = mu*(d-1)/(m*(m-1))`

`a = mu*(m-d)/(m-1)`,

so the valid domain is `mu > 0` and `1 <= d <= m`.

The Fortran regression API accepts an integer response vector and numeric design
matrix, with log or identity link. `m` can be fixed or selected over candidate
orders through `min(max(y),10)`.

## Optimization

The R package switches among `maxLik`, `optim`, and Poisson `glm`. The Fortran
port uses a standalone Nelder-Mead optimizer. For `m>1`, dispersion is
optimized through a logistic transform so trial points stay inside `(1,m)`.
Identity-link fits receive an objective penalty if any fitted mean is
nonpositive.

The final Hessian is calculated numerically in natural `(beta,d)` coordinates
and inverted for a covariance matrix when nonsingular.

Automatic order selection retains the candidate with greatest log-likelihood,
matching the main upstream selection criterion. Formula/model-frame handling
and the intercept-only `ex.solMLE` short-circuit are R orchestration details
and are not separately emulated.

## Likelihood-ratio p-value

For non-Poisson fits, `d=1` is a boundary. The result reports the upstream
half-chi-square(1) tail convention: `0.5 * P(ChiSq_1 >= LR)`.

## RNG

Random generation exactly follows `Poisson(a) + m*Poisson(b)`, but Fortran's
intrinsic RNG supplies the uniform stream, so seeds do not reproduce R's stream
bit-for-bit.
