# API coverage

This port translates numerical algorithms rather than R's formula, S3, zoo,
plotting, printing, and data-frame interfaces. The public Fortran entry point is
the `strucchange` module.

## Directly covered computational areas

| Upstream area | Fortran API | Coverage |
| --- | --- | --- |
| `recresid` | `recursive_residuals` | Recursive residual algorithm and start/end selection |
| Linear regression helpers | `ols_fit`, `inverse_crossprod`, `root_matrix` | SVD least squares, inverse cross-product, symmetric PSD square root |
| `Fstats` | `compute_fstats`, `sup_f_statistic`, `ave_f_statistic`, `exp_f_statistic` | F-statistic sequence and standard aggregations |
| `pvalue.Fstats` | `fstats_pvalue` | Upstream supF/aveF/expF response-surface tables and interpolation |
| `breakpoints` | `compute_breakpoints`, `compute_breakpoint_path`, `best_break_count` | RSS triangle, optimal dynamic programming, RSS/BIC path |
| Segmented regression | `segmented_fit` | Segment coefficients, fitted values, residuals, and segment RSS |
| Breakpoint confidence intervals | `breakpoint_confidence_intervals`, `pargmax_v` | Default upstream asymptotic interval calculation; supports heterogeneous-regressor/error switches |
| OLS/recursive CUSUM | `ols_cusum`, `recursive_cusum` | Empirical fluctuation processes |
| OLS/recursive MOSUM | `ols_mosum`, `recursive_mosum` | Empirical moving-sum processes |
| Recursive/moving estimates | `recursive_estimates_process`, `moving_estimates_process` | Coefficient fluctuation processes with optional rescaling |
| Score processes | `score_cusum`, `score_mosum` | Regression score CUSUM/MOSUM including variance score |
| Generic `gefp` numerical core | `generalized_fluctuation_process` | Cumulative supplied estimating functions, covariance standardization/decorrelation |
| Basic process functionals | `process_max`, `process_range`, `process_max_l2`, `process_mean_l2` | max, range, max-L2, mean-L2 |
| `pvalue.efp` | `efp_pvalue` | Brownian motion, bridge, and increment p-values supported by upstream tables/formulas |
| `supLM` | `sup_lm_statistic`, `sup_lm_pvalue`, `sup_lm_critical_value` | Andrews-style trimmed supLM functional |
| `maxMOSUM` | `max_mosum_statistic`, `max_mosum_pvalue`, `max_mosum_critical_value` | Maximum MOSUM functional |
| `catL2BB` | `cat_l2_bb_statistic`, `cat_l2_bb_pvalue`, `cat_l2_bb_critical_value` | Categorical L2 Brownian-bridge functional |
| Monitoring critical values | `mre_critical_value`, `monitor_me_critical_value`, `monitor_re_critical_value` | Analytical/interpolated critical values |
| Monitoring boundaries | `monitor_ols_cusum_boundary`, `monitor_re_boundary`, `monitor_power_boundary`, `log_plus` | Upstream boundary formulas |

`dp` is re-exported from the public `strucchange` module. It is imported from
`r_kinds` in `rfortran-core`; all maintained real-valued Fortran source uses
`real(dp)` and `_dp` literals.

## Deliberately not translated as public Fortran API

- R formula/model-frame construction, S3 methods, `zoo` time indexes, printing,
  plotting, and `lines` methods.
- R callbacks such as arbitrary `vcov.` functions, `bread`, `estfun`, or custom
  model fit functions. The generalized fluctuation process instead accepts the
  computed score/estimating-function matrix and, optionally, its covariance.
- Stateful R monitoring objects and data-frame update orchestration. Their
  numerical critical values, boundaries, and fluctuation-process building
  blocks are exposed independently.
- `ordL2BB` and `ordwmax`, whose upstream implementation depends on simulated or
  multivariate-normal probability/quantile callbacks from `mvtnorm`. They are
  not silently approximated here.
- Plot-oriented confidence-interval/date formatting and other presentation-only
  helpers.
- Critical-value simulation scripts used upstream to regenerate fixed tables;
  the published tables themselves are retained and tested.

## Notes on breakpoint confidence intervals

The Fortran routine implements the built-in covariance path used when R's
`vcov.` callback is absent. The upstream `het.reg` and `het.err` choices are
available as logical arguments. Arbitrary sandwich/model-specific covariance
callbacks remain outside the package-independent numerical interface.
