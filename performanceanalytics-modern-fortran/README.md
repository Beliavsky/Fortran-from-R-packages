# PerformanceAnalytics Modern Fortran

A modern Fortran translation of the computational parts of the R package `PerformanceAnalytics` 2.1.0. The project uses plain arrays and derived types. Plotting and R object infrastructure are intentionally excluded.

## Implemented numerical areas

### Returns, portfolios, and summaries

- Discrete, logarithmic, difference, cumulative, annualized, excess, relative, centered, and frequency-converted returns
- Geltner unsmoothing and Boudt-style robust cleaning
- Buy-and-hold and periodically rebalanced portfolio accounting
- Beginning/end weights, contributions, turnover, transaction costs, and wealth paths
- Reusable performance-summary statistics and trailing-period outperformance probabilities
- Up, down, and combined capture ratios

### Moments and higher-moment estimation

- Covariance, co-skewness, and co-kurtosis matrices and compact symmetric storage
- Portfolio second, third, and fourth moments and contractions
- EWMA covariance, co-skewness, and co-kurtosis
- Structured M2/M3/M4 targets:
  - Independent unequal marginals
  - Independent equal marginals
  - Observed factor models
  - Latent one-factor/PCA models where applicable
  - Constant-correlation analogues
  - Simaan latent coskewness
  - Central symmetry
- Exact single- and multi-target M2/M3/M4 shrinkage matching the package's analytical loss construction:
  - Exact sample MSE terms
  - Exact target/sample covariance terms
  - Independent, equal-marginal, observed one-factor, and constant-correlation targets
  - Simaan and central-symmetry coskewness targets
  - Multiple observed-factor targets
  - Sixth-order k-statistic correction for unbiased coskewness MSE
  - Returned target vectors, A matrix, b vector, constrained weights, and corrected estimate
- Simpler plug-in shrinkage procedures retained for lightweight use
- Bias-corrected structured independent/equal-marginal coskewness
- Moment Component Analysis for co-skewness and co-kurtosis using higher-order orthogonal iteration
- Nearest Comoment Estimation with PCA initialization, moment-feasibility constraints, optional moment subsets, and identity/diagonal/ridge order weighting

### Tail risk and resampling

- Historical, Gaussian, Cornish-Fisher modified, and lognormal VaR and ES
- Generalized Pareto fitting, VaR, ES, covariance-based uncertainty, and confidence intervals
- Gaussian Monte Carlo VaR and ES for assets and portfolios
- Monte Carlo component VaR and ES
- Kernel portfolio VaR and ES with additive risk contributions
- Moving-block bootstrap standard errors for historical VaR, ES, and Sharpe ratio
- Kupiec and Christoffersen VaR backtests
- Marginal and component risk calculations

### Performance and regression analytics

- Major performance, downside, tail, and drawdown ratios
- CAPM/single-factor regressions and market-timing models
- Rolling and expanding CAPM regressions
- Conditional CAPM with lagged information variables and market interactions
- Rolling and expanding descriptive and risk statistics
- Dynamic correlation arrays

### Applications

```text
analyze_csv data/example_returns.csv 1 2 12
demo_performanceanalytics
portfolio_contributions
advanced_estimators
exact_moment_shrinkage
```

## Build and test

```text
make check
make release
```

`make check` uses GNU Fortran 2018, full runtime checking, floating-point traps, warnings as errors, and the GPL source audit. `make release` repeats all tests and executable workflows at `-O2`.

The project includes `fpm.toml`. `fpm` was unavailable in the validation environment and is therefore not claimed as tested.

## Numerical differences

- Plain arrays replace `xts`, `zoo`, data frames, formulas, and date-aware alignment.
- Periodic rebalancing uses an integer observation interval rather than calendar endpoints.
- NCE uses a deterministic constrained coordinate search. Identity and diagonal/ridge order weighting are available, but the R package's complete `nloptr`/`rgenoud` orchestration and full asymptotic weight-matrix construction are not reproduced exactly.
- `exact_m2_shrinkage`, `exact_m3_shrinkage`, and `exact_m4_shrinkage` reproduce the original analytical target/loss construction. The constrained quadratic program is solved by deterministic projected gradient rather than `quadprog`, so the criterion is preserved but iteration paths may differ.
- The older `shrink_*` and `multi_target_shrink_*` procedures remain lightweight plug-in alternatives.
- Constant-correlation higher-moment targets use symmetry-class averages of standardized sample co-moments.
- GPD estimation uses bounded Nelder-Mead and numerical Hessian/delta uncertainty rather than R's optimizer and profile-root implementation.
- Monte Carlo and bootstrap routines use a reproducible Fortran RNG rather than R's random stream.
- OLS uses the existing small-system pivoted solver in this project.

## Remaining exclusions

- Charts, plotting, colors, legends, GUI elements, and interactive functions
- R formulas, model frames, expression evaluation, S3/S4 classes, attributes, `xts`/`zoo` metadata, and irregular-calendar orchestration
- Formatted report/table presentation; the underlying summary, capture, outperformance, risk, and regression statistics are implemented
- Exact external `RPESE`, `gamlss`, `RobStatTM`, `fit.models`, `quantreg`, `nloptr`, and `rgenoud` behavior
- Exact full NCE optimal/ridge weight matrices and bootstrap ridge-parameter selection
- Packaged datasets and R import helpers

The source audit found analytical correction formulas rather than separate finite-sample lookup-table files for M2/M3/M4 shrinkage. Those formulas are now implemented. Two apparent defects in the original C kernels are corrected and documented: the distinct-index VM3 accumulator now keeps its three mixed moments separate, and the all-distinct constant-correlation CM4 branch no longer divides an already averaged fourth moment by the sample size a second time.

See `API_MAP.md` for the routine-level mapping and `VALIDATION.md` for executed tests.

## License

The original package declares `GPL-2 | GPL-3`, represented here as `GPL-2.0-or-later`. `LICENSE` contains GNU GPL version 2, and every Fortran source, application, example, and test file carries the corresponding SPDX identifier and license notice.
