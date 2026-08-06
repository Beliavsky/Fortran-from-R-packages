# RSDC-fortran

Modern Fortran translation of the computational core of the R package
**RSDC: Regime-Switching Dynamic Correlation Models**.

The library models standardized multivariate observations with regime-specific
correlation matrices and either fixed or covariate-driven transition
probabilities. It is self-contained apart from the vendored, GPL-compatible
`deoptimr-modern-fortran` global optimizer.

## Implemented functionality

- Hamilton filtering and Kim/Hamilton backward smoothing
- Fixed transition matrices and time-varying transition probabilities
  - logistic diagonal persistence for two regimes
  - reference-category softmax for three or more regimes
- Gaussian correlation-model log likelihood and per-observation contributions
- Constant, fixed-transition, and TVTP estimation
- Canonical partial-correlation parameterization, guaranteeing positive-definite
  correlation proposals during global search
- Differential-evolution global search followed by bounded local pattern search
- 70/30 out-of-sample likelihood evaluation
- Data-driven warm-start generation
- Regime ordering by increasing average correlation
- Simulation, Viterbi decoding, filtered/smoothed correlation paths, and
  multi-step-ahead forecasts
- Minimum-variance and maximum-diversification portfolios, with long-only and
  lagged-return options
- Finite-difference scores, Hessian/OPG/sandwich covariance estimates, regime
  diagnostics, parameter-uncertainty correlation bands, and parametric bootstrap
  with parameter draws, covariance, standard errors, and percentile intervals

## Deliberately omitted

Plotting, `ggplot2`, S3/broom methods, R formulas and model frames, printing
infrastructure, datasets, and R parallel-process management are not numerical
algorithms and are omitted. The complete upstream source is retained under
`original/RSDC-master`.

## Build with FPM

```text
fpm build
fpm test
fpm run --example
fpm run
```

The package version is valid SemVer (`0.1.0`), avoiding the invalid
hyphenated-version issue encountered by older FPM releases.

## Build with GNU Fortran

```text
make check
make optimized
make demo
```

On Windows with GNU Fortran, run `run_tests.bat`.

## Minimal filter example

```fortran
use rsdc, only: dp, rsdc_filter_result, rsdc_hamilton

type(rsdc_filter_result) :: out
real(dp) :: y(100, 2), rho(2, 1), p(2, 2)

rho(:, 1) = [0.10_dp, 0.75_dp]
p = reshape([0.95_dp, 0.08_dp, 0.05_dp, 0.92_dp], [2, 2])
call rsdc_hamilton(y, rho, out, pmat=p)
```

Fortran uses `y(time, series)`, `rho(regime, pair)`, and transition matrices
with rows representing the previous state and columns the next state.
Correlation vectors use R's `lower.tri` column-major ordering:
`(2,1), (3,1), ..., (K,1), (3,2), ...`.

## Numerical differences from RSDC

The filter, likelihood, transition maps, forecasts, simulation timing, and
Viterbi recursion follow the upstream formulas. Optimization uses the supplied
DEoptimR-style jDE implementation plus a bounded pattern-search refinement,
rather than R's `DEoptim` plus L-BFGS-B. Numerical derivatives are internal
central finite differences rather than `numDeriv`/`optimHess`. Long-only
portfolio solutions use active-set inverse-covariance formulas rather than
`quadprog`/`Rsolnp`.

See `PORTING_NOTES.md`, `API_MAP.md`, and `VALIDATION.md`.
