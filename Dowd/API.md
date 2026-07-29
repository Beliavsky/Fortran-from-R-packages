# Public API overview

Use the umbrella module:

```fortran
use dowd
```

or import individual modules for tighter namespaces.

## Risk measures

`normal_var`, `normal_es`, `student_t_var`, `student_t_es`, `lognormal_var`,
`lognormal_es`, `log_student_t_var`, `log_student_t_es`, `historical_var`,
`historical_es`, `cornish_fisher_var`, `cornish_fisher_es`, `kernel_var`,
`kernel_es`, `bootstrap_var_es`, `bootstrap_confidence_interval`, `boxcox_var`,
`boxcox_es`, and `spectral_risk_normal`.

Kernel selectors are `kernel_gaussian`, `kernel_box`, `kernel_triangular`, and
`kernel_epanechnikov`.

## Extreme-value methods

`gumbel_var`, `gumbel_es`, `frechet_var`, `frechet_es`,
`generalized_pareto_var`, `generalized_pareto_es`, `hill_estimator`,
`hill_quantile_estimator`, and `pickands_estimator`.

## Portfolio methods

`variance_covariance_var`, `variance_covariance_es`,
`adjusted_variance_covariance_var`, `adjusted_variance_covariance_es`,
`normal_var_hotspots`, `normal_es_hotspots`, `adjusted_var_hotspots`,
`adjusted_es_hotspots`, `pca_prelim`, `pca_var`, and `pca_es`.

## Backtests

`binomial_backtest`, `christoffersen_unconditional_coverage`,
`christoffersen_independence`, `christoffersen_conditional_coverage`,
`lopez_backtest`, `blanco_ihle_backtest`, `jarque_bera_backtest`,
`ks_statistic_normal`, `kuiper_statistic_normal`,
`anderson_darling_statistic_normal`, and `simulate_normal_gof_interval`.

## Options and copulas

Black-Scholes prices, long/short option VaR, simulation ES, American put
binomial/simulation routines, product/Gaussian/Gumbel copulas, sum CDFs, and
copula VaR routines are exposed through `dowd_options` and `dowd_copulas`.

## Utilities

The package exports its internal probability functions, special functions,
sample moments, quantiles, sorting, and random generators because they are
useful for testing and for extending the translated examples.
