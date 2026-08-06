# Source coverage

## `src/stochvolTMB.cpp`

- stationary Gaussian AR(1) latent density: implemented
- Gaussian observation density: implemented
- standardized Student-t observation density: implemented
- standardized skew-normal observation density: implemented
- Gaussian leverage density: implemented
- skew-normal leverage density: implemented
- real-to-(-1,1) parameter transformation: implemented
- Laplace integration formerly supplied by TMB: independently implemented

## `R/simSV.R`

- `sim_sv`: implemented
- `logit`: implemented as `inv_logit_pm1`
- `simulate_parameters`: implemented

## `R/optSV.R`

- objective construction: implemented numerically
- fixed-parameter estimation: implemented
- transformed parameter reporting: exposed through result fields
- latent-mode reporting and standard errors: implemented
- predictive simulation: implemented
- predictive quantile summaries: implemented

## `R/residuals.R`

- deterministic conditional PIT residuals: implemented
- exact TMB `oneStepPredict` correction: not reproduced

## `R/volplot.R`, `R/demo.R`

Omitted because they are plotting or Shiny user-interface code.
