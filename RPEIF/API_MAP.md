# API map

| RPEIF R API | Modern Fortran API | Notes |
|---|---|---|
| `IF` | `influence_series`, `influence_from_data`, `influence_from_nuisance` | String dispatcher with typed results and status codes. |
| `IF.Mean` | `if_mean` or estimator `Mean` | Complete. |
| `IF.SD` | `if_sd` or estimator `SD` | Uses sample SD as the R implementation does. |
| `IF.SemiSD` | `if_semisd` or estimator `SemiSD` | Complete. |
| `IF.VaR` | `if_var` or estimator `VaR` | Gaussian KDE density at empirical quantile. |
| `IF.ES` | `if_es` or estimator `ES` | Complete. |
| `IF.SR` | `if_sr` or estimator `SR` | Source-compatible and corrected modes. |
| `IF.SoR` | `if_sortino` or estimator `SoR` | Constant or mean threshold. |
| `IF.DSR` | `if_downside_sharpe` or estimator `DSR` | Complete. |
| `IF.ESratio` | `if_es_ratio` or estimator `ESratio` | Complete. |
| `IF.VaRratio` | `if_var_ratio` or estimator `VaRratio` | Source-compatible and corrected density modes. |
| `IF.RachevRatio` | `if_rachev_ratio` or estimator `RachevRatio` | Complete. |
| `IF.robMean` | `if_robust_mean` or estimator `robMean` | Uses RobStatTM psi families. |
| `IF.LPM` | `if_lpm` or estimator `LPM` | Orders 1 and 2. |
| `IF.OmegaRatio` | `if_omega_ratio` or estimator `OmegaRatio` | Source-compatible and corrected UPM modes. |
| `nuisParsFn` | `nuisance_parameters_fn` | Returns a typed nuisance-parameter object. |
| `EvaluateShape` | `evaluate_shape` | Returns data only; plotting omitted. |
| `robust.cleaning` | `robust_clean` | Robust winsorization. |
| `LPM` | `lower_partial_moment` | General positive integer order supported by helper. |
| `UPM` | `upper_partial_moment` | Correct and source-compatible sign conventions. |
| `arima(..., order=c(p,0,0))` prewhitening | `ar_prewhiten` | Conditional OLS AR(p), preserving vector length. |
| plotting calls | omitted | Out of scope. |
| `xts` and `zoo` date handling | omitted | Numeric arrays are returned; callers retain dates separately. |
