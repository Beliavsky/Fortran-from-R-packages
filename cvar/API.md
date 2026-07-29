# API reference

All public entities are re-exported by module `cvar`. Floating-point values use `dp = kind(1.0d0)`.

## Risk measures

### Quantile-function input

```fortran
value = var_qf(qf, p_loss [, intercept, slope, transf, status])
value = es_qf(qf, p_loss [, intercept, slope, transf, tol, status])
```

`p_loss` may be a scalar or rank-one array. Array input returns an allocatable array.

### CDF input

```fortran
value = var_cdf(cdf, p_loss [, intercept, slope, transf, tol, status])
value = es_cdf(cdf, p_loss [, intercept, slope, transf, tol, status])
```

The requested quantile is found with an expanding bracket and bisection. `tol` controls quantile or integration accuracy.

### PDF input

```fortran
value = es_pdf(pdf, qf, p_loss [, intercept, slope, transf, tol, status])
```

The quantile function supplies the upper endpoint of the lower-tail partial expectation. The integral from negative infinity is evaluated after a rational variable transformation.

VaR from a PDF alone is not provided, matching the original package, which reports that path as not ready.

### Sample input

```fortran
value = var_sample(sample, p_loss [, intercept, slope, transf, status])
value = es_sample(sample, p_loss [, intercept, slope, transf, status])
```

Supported inputs:

- vector sample and scalar probability: scalar result;
- vector sample and probability vector: vector result;
- matrix sample and scalar probability: one result per column;
- matrix sample and probability vector: matrix with probabilities by columns.

Non-finite sample values are omitted. Quantiles follow R's default type-7 interpolation. Empirical ES averages all observations less than or equal to the interpolated lower-tail cutoff, matching `cvar`.

## Sign convention

The library follows the original package:

```text
VaR_alpha(X) = -q_alpha(X)
ES_alpha(X)  = -E[X | X <= q_alpha(X)]
```

For `Y = intercept + slope * X`, `slope` must be positive. If `transf=.true.`, `Y` is interpreted as a log return and risk is computed for `exp(Y)-1`.

## Distribution utilities

```fortran
normal_pdf(x)
normal_cdf(x)
normal_quantile(p)

student_t_pdf(x, nu)
student_t_cdf(x, nu)
student_t_quantile(p, nu)

std_student_t_pdf(x, nu)
std_student_t_cdf(x, nu)
std_student_t_quantile(p, nu)

ged_pdf(x, nu)
ged_cdf(x, nu)
ged_quantile(p, nu)
ged_scale(nu)
```

The standardized Student-t requires `nu > 2`. The GED uses the `fGarch` variance-one parameterization.

## GARCH(1,1)

### Model construction

```fortran
type(garch11_model) :: model

model = make_garch11(omega, alpha, beta, &
                     innovation=innovation_normal, &
                     shape=5.0_dp, eps0=..., h0=..., eps0sq=...)
```

Innovation constants are:

```fortran
innovation_normal
innovation_std_t
innovation_ged
```

Methods:

```fortran
model%valid()
model%unconditional_variance()
```

A valid model requires `omega > 0`, nonnegative `alpha` and `beta`, and `alpha + beta < 1`.

### Simulation

```fortran
type(garch11_simulation) :: simulation
call simulate_garch11(model, n, simulation [, burnin, seed])
```

Returned fields:

```text
eps     simulated observations
h       conditional variances
eta     standardized innovations
status  status code
```

### Forecasting

```fortran
type(garch11_forecast) :: forecast
call forecast_garch11(model, eps, sigmasq, n_ahead, forecast, &
                      nsim=10000, seed=1234_int64)
```

Returned fields:

```text
eps                  zero conditional-mean forecasts
h                    conditional variance forecasts
plugin_interval      plug-in predictive intervals
simulation_interval  Monte Carlo predictive intervals
simulated_eps         simulated future paths
simulated_h           simulated future variance paths
status                status code
```

The default predictive interval probabilities are 0.025 and 0.975 and can be changed with `lower_probability` and `upper_probability`.

## Status values

```text
cvar_ok
cvar_invalid_probability
cvar_invalid_scale
cvar_invalid_sample
cvar_bracket_failure
cvar_nonconvergence
cvar_invalid_model
cvar_allocation_failure
```

`cvar_status_message(status)` returns a descriptive string.
