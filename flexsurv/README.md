# flexsurv-fortran

Modern free-format Fortran/FPM translation of the computational layer of R
package **flexsurv 2.3.2**.  Plotting, S3 methods, R formula/model-frame
machinery and tidyverse presentation code are intentionally omitted.

The package is standalone: the numerical portions needed from the supplied
`numDeriv`, `deSolve`, and `quadprog` translations are compiled directly from
`src/`, together with the rate-table engine from the translated `relsurv`
package used by advanced `standsurv`. Complete dependency archives are retained
under `vendor/` for provenance.  The supplied survival AFT modules are now compiled for source-style starting values; the complete survival and mstate archives remain under `vendor/` for provenance and integration reference.

## Build

```text
fpm build
fpm test
fpm run --example demo_flexsurv
```

The source uses `dp = kind(1.0d0)` throughout the translated flexsurv modules.

## Implemented computational areas

### Parametric distributions

The native distribution catalogue includes:

- exponential
- Weibull AFT and proportional-hazards parameterizations
- gamma
- log-normal
- Gompertz
- log-logistic
- generalized gamma (Prentice and original parameterizations)
- generalized F (Prentice-style and original parameterizations)

Density, CDF/survival, hazard, cumulative hazard, quantile, random generation,
mean and RMST interfaces are provided where defined.  The generalized gamma,
generalized F, Gompertz and log-logistic formulas follow the upstream C++
implementations.

### `flexsurvreg` numerical core

`fit_flexsurvreg` supports:

- exact events
- right censoring
- left censoring
- interval censoring
- delayed entry / left truncation
- individual right truncation
- observation weights
- background mortality hazards / conditional background survival
- regression effects on any distribution parameter
- fixed transformed parameters
- source-style AFT/moment starting values
- selectable BFGS or Nelder-Mead optimization with optional fallback
- numerical Hessian/covariance calculation through the supplied `numDeriv`
  translation
- AIC, AICc and BIC
- density, survival, hazard, cumulative-hazard and quantile prediction

### User-defined distributions

`custom_distribution` accepts Fortran callbacks for log-density, CDF and
optionally hazard.  `fit_custom_survival` fits the resulting model under the
same censoring/truncation data representation.

### Royston-Parmar splines

Implemented functionality includes:

- restricted cubic spline basis and derivative basis
- `splines2ns`-compatible natural cubic B-spline basis with boundary extrapolation
- hazard, odds and normal scales
- log-time and identity-time scales
- survival, density, hazard, cumulative hazard, quantiles, random generation,
  means and RMST
- spline likelihood fitting
- covariate shifts on the spline predictor
- monotonic initial-value quadratic programming using the supplied `quadprog`
  Goldfarb-Idnani solver

### Mixture models

The mixture layer now includes covariate-dependent multinomial-logit class
membership, known/unknown/partially known event types, arbitrary component
regressions, direct and EM fitting, posterior probabilities, and direct/Louis
covariance estimates.  `fmixmsm` pathway/final probability, mean and quantile
summaries include parameter-bootstrap confidence intervals.

### Fully parametric multi-state models

`flexsurv_multistate` provides:

- time-dependent transition-intensity matrices
- parametric and Royston-Parmar transition models
- transition-probability matrices from the Kolmogorov forward equations
- state-occupation probabilities
- expected length of stay
- path simulation from transition-specific hazards

The ODE calculation uses the supplied `deSolve` Runge-Kutta implementation.
Clock-reset semi-Markov path simulation, predictable time-dependent covariate updates, simulation-based transition matrices, LOS, final-state summaries and parameter-bootstrap intervals are also provided. Shared-regression multi-state models can propagate one joint parameter covariance across all transitions.

### Standardization and utility calculations

Included are weighted standardized survival, hazard, RMST and quantiles,
rate-table-aware all-cause/attributable estimands, contrasts and transformations,
delta/bootstrap uncertainty, Cox-Snell residuals, time-varying hazard ratios,
normal parameter draws, and fractional-polynomial bases/derivatives.

### Right-truncated data

Both computational approaches in the upstream package are represented:

- `survrtrunc_fit`: time-reversed Lynden-Bell/Kaplan-Meier estimator
- `fit_flexsurvrtrunc`: joint and conditional-on-final parametric likelihoods

## v0.3.0 numerical-parity expansion

Version 0.3.0 closes the five numerical targets left after v0.2.0: predictable
time-dependent covariates in semi-Markov simulation, the alternative `splines2ns`
basis, closer flexsurv/survreg starting and optimizer paths, analytic multinomial
Louis derivatives, and cross-transition covariance for a shared multi-state
regression fit. See `docs/TRANSLATION_STATUS.md` for the small remaining
floating-point/interface differences.

## v0.2.0 parity expansion

The second release adds the major numerical workflows that were explicit gaps
in v0.1.0: full covariate-dependent `flexsurvmix` with Louis information,
`fmixmsm` and its parameter-bootstrap summaries, rate-table-aware advanced
`standsurv`, spline coefficient-specific interactions, semi-Markov multi-state
simulation and uncertainty, final-state summaries, AJ comparison helpers, and
fixed-parameter/bootstrap inference for `flexsurvrtrunc`.

The remaining numerical differences are now specialized and are listed in
`docs/TRANSLATION_STATUS.md`.  They chiefly concern predictable time-dependent
covariate updating during semi-Markov simulation, the optional `splines2ns`
basis, and exact R optimizer/initialization micro-behavior.

## License and provenance

The upstream `flexsurv` package is GPL-2-or-later. Because v0.3.0 includes the GPL-3-or-later `splines2ns` compatibility module, the combined v0.3.0 distribution is released under GPL-3.0-or-later; see `LICENSES.md`.  Original upstream source is
preserved under `upstream/flexsurv-2.3.2/`.  See `LICENSES.md` for dependency
licensing and the provenance of translated modules.
