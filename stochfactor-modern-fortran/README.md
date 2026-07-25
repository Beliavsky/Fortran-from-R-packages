# stochfactor-modern-fortran

A combined modern Fortran translation of the computational cores of the R packages
`stochvol` and `factorstochvol`. The factor stochastic-volatility implementation
uses the univariate stochastic-volatility layer in this same project.

This release is a tested numerical translation, not a claim of exact MCMC-chain
identity with the R packages. See `API_MAP.md` and "Explicit limitations" below.

## License

The attached source packages both declare `GPL (>= 2)`. This project therefore
uses `GPL-2.0-or-later` in `fpm.toml` and every Fortran source file. `LICENSE`
contains the GNU General Public License version 2 text.

## Implemented computational features

### Univariate stochastic volatility

- Gaussian stochastic-volatility simulation.
- Variance-standardized Student-t observation errors through inverse-gamma scale
  mixing.
- Leverage through correlation between the observation shock and the next
  log-volatility innovation.
- Linear regression means with a user-supplied design matrix.
- Stationary initialization of the latent log-volatility state.
- Complete-data log likelihood, including the latent-state density, observation
  density, leverage term, Student-t scale density, and supported priors.
- Omori ten-component Gaussian-mixture constants and auxiliary-mixture latent
  state updates for non-leverage models.
- General random-walk Metropolis latent-state updates for leverage models.
- Metropolis updates for `mu`, `phi`, `sigma`, leverage `rho`, Student-t degrees
  of freedom `nu`, and regression coefficients.
- Student-t latent-scale updates with a leverage correction step.
- Gaussian Bayesian regression updates.
- Posterior storage for parameters, regression coefficients, latent states,
  initial latent states, scale mixtures, and acceptance counts.
- Posterior predictive returns, log volatilities, and volatilities.
- Standardized residuals.
- Rolling one-step variance forecasts.
- Log-return conversion with optional demeaning and standardization.

### Factor stochastic volatility

- Simulation from factor-SV models with:
  - User-supplied factor loadings.
  - Idiosyncratic and factor log-volatility processes.
  - Optional homoskedastic components.
  - Optional variance-standardized Student-t factor shocks for simulation.
  - Zero-factor models.
- Dense and sparse default-loading generators and default parameter generators.
- Conditional covariance and correlation paths.
- Individual covariance and correlation elements.
- Ledermann factor bound and exponentially weighted covariance estimation.
- PCA/static-factor initialization.
- Numerical ordering and loading-restriction helpers.
- Gibbs/MCMC estimation with:
  - Conditional Gaussian factor draws.
  - Conditional Gaussian loading draws.
  - Shared univariate SV updates for idiosyncratic and factor volatilities.
  - Lower-triangular identification or a user-supplied logical restriction
    matrix.
  - Sign normalization.
  - Zero-factor support.
- Columnwise Normal-Gamma loading shrinkage with Gamma global updates and GIG
  local-scale updates using a log-scale slice sampler.
- Multi-step draws of log volatilities, covariance matrices, correlation
  matrices, precision matrices, and precision log determinants.
- Woodbury precision and determinant calculations.
- Predictive Gaussian log likelihood using either covariance matrices or the
  Woodbury precision representation.
- Log-mean-exp aggregation over predictive draws.
- Conditional mean and idiosyncratic-volatility prediction.
- Multivariate Normal column densities.
- Sign and factor-order identification.
- Loading-matrix eigenvalue diagnostics.
- Posterior mean covariance and correlation matrices at a selected time.

## Source layout

```text
src/sv_*.f90       shared univariate SV, RNG, statistics, and linear algebra
src/fsv_*.f90      factor-SV types and computations
app/demo_sv.f90    simulated univariate SV fit
app/demo_factor.f90 simulated factor-SV fit
app/fit_csv.f90    CSV command-line application
example/custom_restriction.f90 custom restriction and Normal-Gamma example
test/test_sv.f90   univariate regression tests
test/test_factor.f90 factor-SV regression tests
```

## Build

GNU Fortran with LAPACK and BLAS:

```sh
make debug
make release
```

The debug target uses bounds and runtime checking. Both targets treat warnings as
errors and run all tests and applications.

An `fpm.toml` file is included. `fpm` was not installed in the validation
environment, so the `fpm` build is not claimed as tested.

## CSV application

The input format is:

```text
Date,Series1,Series2,...
2020-01-02,0.001,-0.003,...
```

Fit the first series with univariate SV:

```sh
build/fit_csv data/example.csv sv
```

Fit a factor-SV model with two factors:

```sh
build/fit_csv data/example.csv factor 2
```

The application deliberately exposes a small reproducible interface. Library
users should configure `sv_mcmc_options`, `fsv_options`, priors, restrictions,
and prediction settings directly.

## Explicit limitations

The following are not claimed:

- Exact reproduction of R random-number streams or chain-by-chain posterior
  draws.
- The original optimized C++ fast sampler, adaptive MALA/RWMH machinery, ASIS,
  parameter expansion, or every `expert` tuning option from `stochvol`.
- Every R prior class and every prior combination. The Fortran `sv_prior` type
  implements a compact supported set of Normal, Beta, and Gamma/exponential-like
  prior terms.
- Model-misspecification reweighting and all specialized summary/runtime fields.
- Multiple-chain orchestration, `coda` objects, parallel execution, formula
  interfaces, dates, `ts` metadata, S3 methods, printing, summaries, and plots.
- The exact shallow/deep interweaving sampler used by `factorstochvol`.
- Exact `factanal` initialization. `static_factor_initialize` is a tested
  PCA-based numerical analogue.
- Exact `GIGrvg` random draws. Local Normal-Gamma scales use a tested log-scale
  slice sampler targeting the GIG density.
- Student-t factor innovations in factor-model estimation. They are available
  only in factor-model simulation in this release.
- Complete equivalence for all original validation, default-selection, and
  storage options.

These omissions are computationally important. This release should be described
as a broad combined numerical translation, not complete option-for-option parity.
