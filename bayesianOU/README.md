# bayesianOU-fortran

Modern Fortran/FPM translation of the computational algorithms in
`bayesianOU` 0.2.0, a package for nonlinear Ornstein-Uhlenbeck models with
cubic drift, stochastic volatility, Student-t innovations, nested latent
production prices, model diagnostics, multiple-imputation pooling, and a
separate geometry-adaptive HMC engine.

## Included numerical functionality

- Single-level nonlinear OU/TMG model.
- Two-level market-price to latent-production-price model.
- Three-level extension with a value anchor.
- Cubic restoring drift, COM mean effect, TMG-modulated reversion, Student-t
  innovations, stochastic-volatility state filtering, and level-richness
  switches.
- Training-only standardization, weighted COM statistics, first-principal-
  component common factors, and TMG orthogonalization.
- Native initialization by nonlinear least squares, OLS, and conditional SV
  filtering.
- Empirical-Bayes COM hierarchy shrinkage.
- Blocked random-walk Metropolis posterior sampling of structural parameters.
- Pointwise log likelihoods and a self-contained PSIS-style LOO calculation.
- Single-level and nested conditional out-of-sample forecast metrics.
- R-hat, effective sample size, divergence counts, and kappa stability evidence.
- Rubin multiple-imputation pooling.
- Beta(TMG), drift-decomposition, SV-scale, latent-mean, and accounting outputs
  that replace the numerical content underlying the R plots.
- Generic Euclidean and SoftAbs Riemannian static HMC, including generalized
  implicit leapfrog integration and E-BFMI.

The project has no Stan, R, or C++ runtime dependency. BLAS and LAPACK are used
for dense linear algebra.

## Important backend difference

The original package delegates the OU posterior to Stan/NUTS and jointly
samples all latent stochastic-volatility and latent-production states. This
Fortran port is self-contained and therefore uses a different inference
backend:

1. nonlinear/linear conditional initialization;
2. a robust AR(1) log-squared-residual SV filter;
3. anchor-conditioned latent production paths;
4. empirical-Bayes hierarchy shrinkage; and
5. blocked random-walk Metropolis sampling of structural parameters.

The model equations, transformations, priors, likelihood branches, generated
pointwise log likelihoods, and downstream calculations are translated, but the
posterior draws are not bitwise or algorithmically identical to Stan/NUTS.
Use the separate `ou_geom_hmc` API when an exact user-supplied differentiable
target and static HMC are desired. See `PORTING.md`.

## Build

With FPM:

```text
fpm build
fpm test
fpm run
fpm run --example nested_three_level
fpm run --example geometry_hmc
```

GNU Fortran without FPM:

```text
./build_gfortran.sh strict
./build_gfortran.sh release
```

The GNU script links LAPACK and BLAS.

## Minimal use

```fortran
use bayesianou

type(ou_input)      :: data
type(ou_options)    :: options
type(ou_fit_result) :: fit

! Allocate and fill data%y, data%x, data%tmg, data%com, data%capital.
options%n_levels = 1
options%chains = 4
options%iterations = 2000
options%warmup = 1000

call fit_ou_nonlinear_tmg(data, options, fit)
print *, fit%summary%kappa
print *, fit%diagnostics%loo%elpd_loo
```

For `n_levels >= 2`, also allocate `data%gprime`. For `n_levels == 3`, allocate
`data%value_anchor`.

## Project layout

- `src/`: library modules.
- `app/`: runnable single-level demonstration.
- `example/`: nested and geometry examples.
- `test/`: strict numerical and end-to-end tests.
- `original/`: unmodified original R/Stan package tree.
- `API.md`: public Fortran API mapping.
- `PORTING.md`: design choices and differences.
- `TESTING.md`: validation commands and coverage.
- `REFERENCE_GENERATION.md`: independent fixed references.

## License

MIT. See `LICENSE` and `NOTICE.md`.
