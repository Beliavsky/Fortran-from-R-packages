# SVDNF for Modern Fortran

A dependency-free modern Fortran/FPM translation of the numerical algorithms
in the R package **SVDNF 0.1.11**, which implements Kitagawa discrete nonlinear
filtering for stochastic-volatility models.

## Implemented numerical scope

- Kitagawa discrete nonlinear filtering and log-likelihood evaluation.
- Heston, Bates, Duffie-Pan-Singleton, Taylor, Taylor with leverage,
  Pitt-Malik-Doucet, and stochastic-volatility CAPM dynamics.
- User-defined state drift and diffusion callbacks.
- User-defined jump-count probability and simulation callbacks.
- Automatic volatility, jump-count, and jump-size grids.
- Poisson/Bernoulli jump counts, gamma volatility jumps, and normal return jumps.
- State and return simulation with optional deterministic seeds.
- Filtering and one-step prediction distributions.
- Filtering/prediction percentiles.
- Factor-adjusted returns.
- Native Nelder-Mead maximum-likelihood estimation for built-in models.
- Callback-based parameter mapping for custom-model estimation.
- Numerical Hessians, inverse-Hessian standard errors, and convergence metadata.
- Monte Carlo forecasts for latent volatility and returns.
- Native normal, gamma, Poisson, and binomial probability support.

The R S3 printing/plotting layer, `xts`/`zoo` date adapters, and graphics are not
compiled. Numerical arrays used by those methods are returned directly.

## Build

```text
fpm build
fpm test
fpm run svdnf_demo
fpm run --example built_in_models
fpm run --example filter_and_forecast
fpm run --example estimate_taylor
```

GNU Fortran validation scripts are also included:

```text
sh scripts/validate.sh
```

On Windows:

```bat
scripts\validate.bat
```

## Basic filtering example

```fortran
use svdnf

type(svm_dynamics) :: dynamics
type(simulation_result) :: simulated
type(filter_result) :: filtered

dynamics = dynamics_svm('Heston', mu=0.04_dp, kappa=3.0_dp, &
  theta=0.04_dp, sigma=0.35_dp, rho=-0.65_dp)

simulated = model_simulate(dynamics, 250, initial_volatility=0.04_dp, &
  seed=12345)

filtered = dnf_filter(dynamics, simulated%returns, n=40)

if (.not. filtered%ok) error stop trim(filtered%message)
print *, filtered%log_likelihood
```

## Built-in parameter order

The maximum-likelihood interfaces use these parameter vectors:

| Model | Parameter order |
|---|---|
| Heston | `mu, kappa, theta, sigma, rho` |
| Bates | `mu, alpha, delta, omega, kappa, theta, sigma, rho` |
| Duffie-Pan-Singleton | `mu, alpha, delta, rho_z, nu, omega, kappa, theta, sigma, rho` |
| Taylor | `phi, theta, sigma` |
| Taylor with leverage | `phi, theta, sigma, rho` |
| Pitt-Malik-Doucet | `phi, theta, sigma, rho, delta, alpha, p` |
| CAPM-SV | factor coefficients, then `phi, theta, sigma` |

Use `parameter_vector`, `set_parameter_vector`, and `model_parameter_names` to
avoid hard-coding these layouts.

## Custom dynamics

A custom model supplies four callbacks:

```fortran
call set_custom_dynamics(dynamics, mu_y, sigma_y, mu_x, sigma_x, &
  mu_y_parameters=..., sigma_y_parameters=..., &
  mu_x_parameters=..., sigma_x_parameters=...)
```

Custom maximum likelihood is available through `dnf_optimize_custom`, which
accepts a typed parameter-setter callback. This replaces R's run-time function
argument introspection with compile-time interfaces.

## Main result types

- `svm_dynamics`
- `grid_type`
- `simulation_result`
- `filter_result`
- `percentile_result`
- `forecast_result`
- `optimization_result`

All public numerical arrays use `real(dp)`, where `dp = kind(1.0d0)`.

## Documentation

- `COVERAGE.md`: mapping from upstream routines to Fortran APIs.
- `PORTING_NOTES.md`: design differences and corrected behavior.
- `VALIDATION.md`: compiler flags and validation results.

## License

GPL-3.0-only. See `LICENSE`, `NOTICE`, and the retained upstream source.
