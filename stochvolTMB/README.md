# stochvolTMB-fortran

Modern Fortran 2018 translation of the computational code in the R package
`stochvolTMB` 0.3.0.

The library fits stochastic-volatility models by minimizing a Laplace
approximation to the likelihood after integrating out the latent AR(1)
log-volatility process. It is self-contained and does not require R, TMB,
RcppEigen, BLAS, or LAPACK.

## Implemented models

- Gaussian stochastic volatility
- Standardized Student-t stochastic volatility
- Standardized skew-normal stochastic volatility
- Gaussian stochastic volatility with leverage
- Skew-normal stochastic volatility with leverage

## Main API

```fortran
use stochvoltmb

type(sv_rng_state) :: rng
type(sv_parameters) :: params
type(sv_simulation) :: sim
type(sv_fit_result) :: fit
type(sv_prediction) :: pred

params%sigma_y = 0.3_dp
params%sigma_h = 0.25_dp
params%phi = 0.92_dp

call rng%seed(12345_8)
call sim_sv(params, 500, sv_gaussian, rng, sim)
call estimate_parameters(sim%y, sv_gaussian, fit)
call predict_sv(fit, 10, 5000, rng, pred)
```

Important public procedures include:

- `sim_sv`
- `estimate_parameters`
- `get_nll` and `laplace_nll`
- `latent_mode` and `joint_nll`
- `simulate_parameters`
- `predict_sv`
- `summarize_prediction`
- `standardized_residuals` and `one_step_residuals`
- parameter transformations and model-name conversions

## Building

With FPM:

```text
fpm test
fpm run --target demo_stochvoltmb
```

With GNU Make:

```text
make checked
make optimized
make MODE=optimized demo
```

The checked build enables bounds checking, floating-point traps, backtraces,
and warnings as errors.

## Numerical implementation

TMB is not translated or linked. The replacement engine uses:

- sparse Newton iterations for the latent mode;
- a symmetric tridiagonal latent Hessian;
- an O(n) tridiagonal solve and log determinant;
- the standard Laplace correction;
- Nelder-Mead optimization of transformed fixed parameters;
- finite-difference observed-information matrices;
- exact diagonal extraction from the inverse latent tridiagonal Hessian.

Gaussian and Student-t latent derivatives are analytic. Skew-normal and
leverage observation derivatives are evaluated locally by centered finite
differences, preserving the O(n) cost of each latent Newton step.

## Source compatibility notes

The model equations and parameter transformations follow the upstream C++ and
R code. The leverage likelihood intentionally omits the final observation,
because the upstream TMB objective needs `h(i+1)` and evaluates leverage
observations only for `i < N-1` in zero-based indexing.

The R documentation for the Student-t simulator mentions a scale involving 2,
but the executable R code and TMB likelihood use `sqrt((df-2)/df)`. This port
follows the executable code.

`one_step_residuals` provides conditional probability-integral-transform
residuals evaluated at the fitted latent mode. It is a deterministic diagnostic,
not a bit-for-bit reimplementation of TMB's computationally intensive
`oneStepPredict` correction.

Parameter and latent standard errors are numerical Laplace approximations.
They will not be bit-for-bit identical to TMB automatic-differentiation output.

## Omitted R-only features

- ggplot2 volatility and forecast plotting
- Shiny demonstration application
- S3 classes, print methods, and data.table formatting
- R formula/model-frame infrastructure
- TMB automatic differentiation and `sdreport` objects

## License

GPL-3.0-only, matching the upstream package. The original package source and
archive are retained under `upstream/`. The skew-normal formulas and generator
reuse GPL-3-compatible work from the previously completed `sn` Fortran port;
that archive is also retained for provenance.
