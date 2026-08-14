# KrigInv R-to-Fortran API mapping

| R / dependency operation | Fortran counterpart | Notes |
|---|---|---|
| `DiceKriging::km` | `fit_krig_model` | Native MLE/PMLE/LOO GP fitting |
| fixed `km` covariance parameters | `init_dice_krig_model` | Uses DiceKriging covariance conventions |
| v0.1 fixed model | `init_krig_model` | Backward-compatible legacy mode |
| `update(..., cov.reestim=...)` | `update_krig_model(..., cov_reestimate=...)` | Genuine covariance/variance refit in Dice-backed mode |
| `model@param.estim` | `model%param_estim` | Drives default EGI covariance re-estimation |
| `bichon_optim` | `bichon_optim` | Native formula |
| `computeQuickKrigcov` | `compute_quick_krigcov` | Posterior cross-covariance |
| `computeRealVolumeConstant` | `compute_real_volume_constant` | Native bivariate-normal integration |
| `EGI` | `egi` | Sequential workflow; supports covariance/trend/nugget re-estimation controls |
| `EGIparallel` | `egi_parallel` | Batch workflow; no R cluster layer |
| `excursion_probability` | `excursion_probability` | One or multiple thresholds |
| `integration_design` | `integration_design` | Sobol, MC, SUR, Vorob, JN, TIMSE, IMSE |
| `jn_optim_parallel` | `jn_optim_parallel` | Native batch criterion |
| `jn_optim_parallel2` | `jn_optim_parallel2` | Candidate-plus-existing-batch wrapper |
| `max_infill_criterion` | `max_infill_criterion` | Discrete or differential-evolution optimization |
| `max_futureVol_parallel` | `max_futurevol_parallel` | Includes conservative-alpha fallback |
| `max_sur_parallel` | `max_sur_parallel` | Also handles JN real-volume-variance mode |
| `max_timse_parallel` | `max_timse_parallel` | TIMSE and IMSE |
| `max_vorob_parallel` | `max_vorob_parallel` | Vorob'ev / conservative alpha |
| `precomputeUpdateData` | `precompute_update_data` | Compatibility object |
| `predict_nobias_km` | `predict_nobias_km` | Native SK/UK prediction; Dice backend when fitted |
| `predict_update_km_parallel` | `predict_update_km_parallel` | Batch conditional update |
| `ranjan_optim` | `ranjan_optim` | Native formula |
| `sur_optim_parallel` | `sur_optim_parallel` | Native criterion |
| `sur_optim_parallel2` | `sur_optim_parallel2` | Wrapper |
| `timse_optim_parallel` | `timse_optim_parallel` | Native criterion |
| `timse_optim_parallel2` | `timse_optim_parallel2` | Wrapper |
| `tmse_optim` | `tmse_optim` | Native formula |
| `tsee_optim` | `tsee_optim` | Native formula |
| `vorob_optim_parallel` | `vorob_optim_parallel` | Native criterion |
| `vorob_optim_parallel2` | `vorob_optim_parallel2` | Wrapper |
| `vorob_threshold` | `vorob_threshold` | Boundary-safe interpolation |
| `vorobVol_optim_parallel` | `vorobvol_optim_parallel` | Native future-volume criterion |
| `vorobVol_optim_parallel2` | `vorobvol_optim_parallel2` | Wrapper |
| `print_uncertainty*` | omitted | Plotting/graphics |

## DiceKriging fit controls exposed by KrigInv-fortran

`km_control` is re-exported by the umbrella `kriginv` module.  It carries multistart count, candidate population size, iteration limit, tolerance, gradient use, and the nugget-mixture upper bound used by the translated DiceKriging fitter.

`fit_krig_model` also accepts fixed/initial covariance parameters, process variance, trend coefficients, nugget/nugget estimation, heteroskedastic noise, isotropic mode, explicit bounds, initial optimizer parameters, scaling axes, estimation method, and SCAD penalty parameter.
