# MSGARCH Modern Fortran

A modern Fortran computational translation of the `MSGARCH` R package. The project provides Markov-switching and static-mixture conditional-volatility models without R classes, plotting, `zoo` indexes, or other R-specific infrastructure.

The original package metadata is retained under `reference/`. The original package declares `GPL (>= 2)`, preserved here as `GPL-2.0-or-later` in `LICENSE`, `fpm.toml`, and every Fortran source file.

## Implemented model families

Each regime may use one of:

- `sARCH`
- `sGARCH`
- `gjrGARCH`
- `eGARCH`
- `tGARCH`

Each regime may use one of six standardized innovation laws:

- `norm`
- `std`
- `ged`
- `snorm`
- `sstd`
- `sged`

A specification may contain one regime, a Markov transition matrix, or static mixture weights.

## Implemented computational features

- Specification construction and validation
- Parameter packing, unpacking, bounds, bounded maps, and simplex maps
- Fixed parameters and regime-constant parameter groups during ML estimation
- Variance-targeting intercept formulas used by the original starting-value calculations
- Regime recursion, simulation, likelihood evaluation, and unconditional regime variances
- Hamilton predicted and filtered probabilities
- Kim-style smoothed probabilities
- Viterbi state decoding
- Transition-matrix powers and future state probabilities
- In-sample conditional volatility, density, CDF, and PIT
- One-step analytical predictive density and CDF
- Multi-step predictive density by Gaussian kernel estimation and empirical predictive CDFs
- Monte Carlo mean and volatility forecasts
- Simulation-based unconditional volatility
- Conditional and unconditional simulation
- VaR and Expected Shortfall, in sample and out of sample
- Cumulative-horizon risk calculations
- Bounded maximum-likelihood estimation
- Numerical Hessian, covariance matrix, standard errors, AIC, and BIC
- Adaptive random-walk Metropolis sampling with configurable Gaussian priors
- Burn-in, thinning, posterior means, posterior standard deviations, and DIC
- Optional regime relabeling by unconditional variance for homogeneous specifications
- Posterior state probabilities, volatility, PIT, predictive density/CDF, unconditional volatility, and VaR/ES
- Gaussian HMM and Gaussian mixture EM routines used as numerical initialization/diagnostic utilities

## Building and testing

GNU Fortran:

```text
make check
make release-check
```

or directly:

```text
./run_checks.sh debug
./run_checks.sh release
```

On Windows with GNU Fortran:

```text
run_checks.bat debug
run_checks.bat release
```

The debug configuration enables full GNU Fortran runtime checking. Both configurations treat warnings as errors.

An `fpm.toml` manifest is included. It was not validated because `fpm` was unavailable in the translation environment.

## Applications

Run the simulation/filtering/ML demonstration:

```text
build/debug/demo_msgarch
```

Run the MCMC example:

```text
build/debug/mcmc_example
```

Fit the final numeric column of a CSV file:

```text
build/debug/fit_csv data/returns.csv single sGARCH norm
build/debug/fit_csv data/returns.csv markov sGARCH norm gjrGARCH std
build/debug/fit_csv data/returns.csv mixture sGARCH norm eGARCH ged
```

The CSV reader accepts a single numeric column or a header/date column followed by the numeric series in the final field.

## Minimal library example

```fortran
program example
   use msgarch
   implicit none
   type(msgarch_spec) :: spec
   type(simulation_result) :: simulation
   type(filter_result) :: filtered

   call seed_rng(12345)
   spec = create_spec([character(len=12) :: 'sGARCH', 'gjrGARCH'], &
      [character(len=8) :: 'norm', 'std'])
   simulation = simulate_msgarch(spec, 1000, 1)
   filtered = hamilton_filter(spec, simulation%draw(1,:))
   write(*,'(a,f12.4)') 'log likelihood: ', filtered%loglik
end program example
```

## Numerical differences from the R package

This is a computational translation, not an ABI or byte-for-byte port.

- Bounded Nelder-Mead replaces the package's R optimization orchestration.
- The MCMC routine is an adaptive independent-coordinate random-walk Metropolis analogue, not the exact Vihola adaptive C++ sampler.
- Gaussian priors are directly configurable; the complete R prior-control object system is not reproduced.
- Multi-step predictive densities use a tested Gaussian KDE with a Silverman bandwidth rather than exact reproduction of R's `density()` bandwidth selection.
- Gaussian HMM and mixture EM utilities are deterministic numerical analogues of the package's starting-value machinery.
- Random-number streams do not reproduce R or Rcpp streams.

These differences can produce different fitted endpoints and posterior chains even when the model equations are the same.

## Explicit exclusions

- S3 classes and methods
- `coda` objects and diagnostics
- R formulas and list-based control objects
- `zoo` and `ts` metadata
- Plotting, fan charts, printing, and summaries
- Parallel chains
- Arbitrary R callback samplers
- Exact Vihola/adaptMCMC covariance adaptation
- Complete original automatic starting-value search orchestration
- Exact R optimizer, integration, KDE, or RNG equivalence

See `API_MAP.md` and `VALIDATION.md` for the exact translated surface and test coverage.
