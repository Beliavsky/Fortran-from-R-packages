# Translation status

Source: R package `relsurv` 2.3-3.

## Implemented computational areas

### Population mortality and rate tables

- `pystep.c` / `pystep2.c` rate-table stepping
- `exps.c` expected population survival
- `expprep2`-style expected-survival and summary facades
- population cumulative-hazard increments
- generic rate-table construction and joining
- `transrate` annual-survival/hazard conversion
- `transrate.hld` Human Life-Table parser and conversion
- `transrate.hmd` Human Mortality Database parser and conversion

### Non-parametric relative survival

- Pohar-Perme net survival
- Ederer I and II
- Hakulinen
- Kaplan-Meier and Fleming-Harrington cumulative transforms
- standard-error/confidence-limit calculations in `rs_surv_result`
- `netwei` and `netfastpinter2` numerical accumulators

### Comparison and planning

- `cmp.rel` disease/population crude mortality and variance recursions
- `rs.diff` grouped net-survival score test and covariance
- `nessie` expected numbers alive and expected survival-time integration

### Regression models

- Aalen additive-hazard kernels translated from `aalen_beta.cpp`
- relative/additive Aalen increments
- direct and piecewise additive relative-survival likelihood fitting
- `rsadd` EM branch for unknown cause of death
- Ederer-II-based automatic EM bandwidth selection
- smoothed/unsmoothed EM excess baseline hazards
- grouped binomial and Poisson `rsadd` GLM likelihoods
- transformed-time Cox fitting corresponding to `rstrans`
- multiplicative Cox fitting with population-hazard offset corresponding to `rsmul`
- Cox fitting uses the supplied Fortran `survival` implementation

### Residuals and proportionality diagnostics

- `residuals.rsadd` Schoenfeld/partial-residual numerical construction
- event-specific and summed residual covariance matrices
- `rs.br` Brownian-bridge diagnostics
- maximum-deviation and Cramer-von Mises tests
- additive/rsadd and Cox tied-event conventions
- `rs.zph` identity, rank, log, and KM time transforms
- per-event (`each`) and summed covariance scaling

### Smoothing

- ordinary Epanechnikov smoother
- exact `epa` adaptive-bandwidth segmentation
- asymmetric boundary kernels
- left-predictable Epanechnikov kernel

### Years-lost calculations

- default observed-minus-population years difference
- YL2013 excess-failure integration
- YL2017 observed/population integration
- Greenwood probability standard errors
- upstream Greenwood area variance
- bootstrap column variance for probabilities, areas, and years estimates
- normal confidence intervals

### Other utilities

- counting-process follow-up splitting (`survsplit` numerical core)
- monotone inverse-time interpolation (`invtime` numerical core)
- competing-risk cumulative-incidence helper

## Intentionally modernized interfaces

- R formulas, model frames, `Surv`, `ratetable`, factors, and S3 objects are
  replaced by explicit typed arrays and result types.
- `expprep2` takes a `ratetable_type` explicitly rather than an R formula and
  model frame.
- `rs.br` and `rs.zph` expose their numerical diagnostic results but omit their
  graphical plotting interfaces.
- `years` accepts bootstrap replicate curves explicitly. The numerical
  resampling aggregation and variance formulas are translated; R's report/data-
  frame/plot orchestration is not reproduced.
- `transrate.hld` / `transrate.hmd` return native `ratetable_type` objects rather
  than R ratetable attributes/factors.

## Remaining differences

The important statistical algorithms targeted for v0.2.0 are implemented.
Remaining differences are primarily R runtime/interface semantics:

- exact R formula/model-frame and S3 `predict`, `residuals`, `survfit`, print, and
  summary dispatch
- graphical routines (`plot.*`, `plot_f`, `plot_years`, bandwidth plots)
- exact R `NA`, recycling, factor labels, date classes, and warning behavior
- automatic subject-resampling/report construction around `years`; native code
  instead takes the resulting bootstrap curves explicitly
- byte-for-byte parser behavior for every malformed or nonstandard HLD/HMD file;
  the documented standard formats and numerical conversions are implemented

The complete upstream source remains under `upstream/relsurv` for provenance.
