# Porting notes

## Source package

- Package: `cvar`
- Version: `0.6`
- Original author: Georgi N. Boshnakov
- Original license: `GPL (>= 2)`
- Port license: `GPL-2.0-or-later`

## API mapping

| R entry point | Fortran equivalent |
|---|---|
| `VaR(..., dist.type="qf")` | `var_qf` |
| `VaR_qf` | `var_qf` |
| `VaR(..., dist.type="cdf")` | `var_cdf` |
| `VaR_cdf` | `var_cdf` |
| `VaR.numeric` | `var_sample` with vector input |
| `VaR.matrix` | `var_sample` with matrix input |
| `ES(..., dist.type="qf")` | `es_qf` |
| `ES(..., dist.type="cdf")` | `es_cdf` |
| `ES(..., dist.type="pdf")` | `es_pdf` |
| `ES.numeric` | `es_sample` with vector input |
| `ES.matrix` | `es_sample` with matrix input |
| `GarchModel` | `make_garch11` and `garch11_model` |
| `sim_garch1c1` | `simulate_garch11` |
| `predict.garch1c1` | `forecast_garch11` |

## R facilities replaced

- S3 classes and method dispatch were replaced by generic Fortran interfaces and derived types.
- R's `...` parameter recycling was replaced by explicit callbacks and typed arguments.
- R function names and unevaluated calls were replaced by procedure arguments.
- R's `integrate()` was replaced by adaptive 15-point Gauss-Kronrod integration.
- `gbutils::cdf2quantile` was replaced by an internal expanding-bracket bisection solver.
- R's default `quantile()` behavior was reproduced with type-7 interpolation.
- R RNG-state save/restore semantics were replaced by an explicit deterministic `rng_state`.
- Conditional loading of `fGarch` was replaced by self-contained standardized-t and GED implementations.

## Intentional numerical differences

### Random-number streams

The Fortran port uses a Park-Miller generator plus Box-Muller, Marsaglia-Tsang gamma sampling, and derived Student-t/GED generators. A numeric seed is reproducible within this library, but it does not reproduce R's Mersenne-Twister stream or the saved `.RDS` test fixtures.

### CDF inversion

The original package delegates CDF inversion to `gbutils::cdf2quantile`. The port expands a bracket around zero and then bisects. Results normally agree near machine precision for smooth distributions, but the exact last bits and failure behavior can differ.

### Infinite-tail integration

For quantile/CDF ES, the port applies `p = p_loss * u^2` before adaptive integration, smoothing the integrable endpoint singularity at probability zero. PDF ES maps the negative-infinite interval to `(0,1)`. These transformations improve robustness without changing the mathematical quantity.

### Corrected transformed sample ES

The R implementation computes an intercept/slope-adjusted VaR cutoff but compares it directly with untransformed sample observations. That is inconsistent when the location or scale differs from `(0,1)`. The port selects the lower tail in the original sample scale and then transforms each selected observation.

### Corrected transformed PDF ES location sign

The R `dist.type="pdf", transf=TRUE` branch uses `exp(-intercept + slope*x)`, while the quantile/CDF branches and the stated transformation require `exp(intercept + slope*x)`. The port uses the mathematically consistent positive intercept sign.

## GARCH conventions

The recursion is preserved:

```text
h_t   = omega + alpha * eps_(t-1)^2 + beta * h_(t-1)
eps_t = sqrt(h_t) * eta_t
```

If initial values are absent, both the prior squared observation and prior variance use `omega / (1-alpha-beta)`, as in the R source. Forecast variances use

```text
h_(T+1|T) = omega + alpha*eps_T^2 + beta*h_T
h_(T+k|T) = omega + (alpha+beta)*h_(T+k-1|T), k > 1.
```

## Out of scope

The port does not reproduce R call objects, environments, expression mutation, deprecated argument warnings, package namespace behavior, plotting examples, or external `fGarch` model fitting. The package itself does not fit GARCH parameters; it simulates and forecasts from supplied parameters, and the Fortran port has the same scope.
