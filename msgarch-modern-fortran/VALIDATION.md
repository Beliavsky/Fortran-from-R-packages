# Validation

## Environment

- GNU Fortran 14.2.0
- GNU/Linux x86-64
- Fortran 2018 mode

No external numerical library is required.

## Compiler configurations

Debug:

```text
-std=f2018 -ffree-line-length-none -O0 -g
-Wall -Wextra -Wimplicit-interface -Werror
-fcheck=all -fbacktrace
```

Release:

```text
-std=f2018 -ffree-line-length-none -O2
-Wall -Wextra -Wimplicit-interface -Werror
-fbacktrace
```

## Commands

```text
./test/check_license.sh
./run_checks.sh debug
./run_checks.sh release
```

## Test groups

### Distribution and regime-model tests

- Density integration for all six innovations
- CDF/quantile inversion
- Simulated standardized mean and variance
- Simulation/filtering for all five volatility recursions
- Parameter counts, bounds, packing, and unpacking
- Original variance-targeting formulas
- Markov and mixture specification validation
- Bounded and simplex parameter maps
- Homogeneous-regime detection and unconditional-variance relabeling
- Transition-matrix permutation during relabeling

### Filter, simulation, forecast, PIT, and risk tests

- Hamilton predicted and filtered probabilities
- Smoothed probabilities
- Viterbi paths
- One-step analytical PDF and CDF
- In-sample PDF, CDF, PIT, and conditional volatility
- Conditional simulation
- Mean and volatility forecasts
- Multi-step KDE density and empirical CDF forecasts
- Simulation-based unconditional volatility
- Transition powers and state forecasts
- VaR and Expected Shortfall, including retained paths
- In-sample one-step risk calculations

### Estimation, MCMC, posterior, and HMM tests

- ML likelihood improvement from a perturbed start
- Numerical Hessian and covariance output
- All-fixed parameter evaluation
- Regime-constant tied parameters
- Single-regime extraction
- MCMC burn-in, thinning, acceptance, and finite draws
- Posterior means and DIC
- Posterior state probabilities
- Posterior predictive PDF/CDF, volatility, PIT, unconditional volatility, and risk
- Gaussian HMM EM, transition probabilities, and Viterbi path
- Gaussian mixture EM

## Application checks

- `demo_msgarch`
- `mcmc_example`
- `fit_csv` in single-regime mode
- `fit_csv` in Markov-switching mode
- `fit_csv` in static-mixture mode

Application execution is checked. The sample Markov and mixture fits are intentionally small command-line smoke tests and are not used as convergence benchmarks.

## Interpretation

Passing tests establish that the implemented Fortran procedures compile, execute, satisfy the stated numerical identities and invariants, and work together through the tested workflows. They do not establish exact equality with R optimizer endpoints, posterior chains, KDE values, or random-number streams.

## Linker diagnostic

GNU ld reports that `msgarch_estimation.o` requires an executable stack because internal Fortran objective procedures are passed as optimizer and Hessian callbacks. Linking and all executable tests succeed. Compiler diagnostics remain treated as errors; this is a linker diagnostic and is recorded here rather than hidden.
