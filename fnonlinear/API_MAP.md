# API map

| Original routine | Modern Fortran procedure | Status |
|---|---|---|
| `tentSim` | `tent_sim` | Implemented and tested |
| `henonSim` | `henon_sim` | Implemented and tested |
| `ikedaSim` | `ikeda_sim` | Implemented and tested |
| `logisticSim` | `logistic_sim` | Implemented and tested |
| `lorentzSim` | `lorentz_sim` / `lorenz_sim` | Implemented and tested |
| `roesslerSim` | `roessler_sim` / `rossler_sim` | Implemented and tested |
| `.rk4` | `rk4_integrate_times` | Implemented and tested |
| `mutualPlot` numerical output | `mutual_information_curve` | Implemented and tested |
| `.embeddPSR` | `delay_embed` | Implemented and tested |
| `falsennPlot` numerical output | `false_nearest_neighbors` | Implemented and tested |
| `recurrencePlot` numerical relation | `recurrence_matrix` | Implemented and tested |
| `separationPlot` numerical output | `space_time_separation` | Implemented and tested |
| `lyapunovPlot` numerical output | `lyapunov_stretching` | Implemented and tested |
| `.find.nearest` | `find_k_nearests` | Implemented and tested |
| `.follow.points` | internal/public neighbor-following engine | Implemented and tested through Lyapunov paths |
| `.lyapunovFit` | `lyapunov_linear_fit` | Implemented and tested |
| `.C2` | `correlation_integral` | Implemented and tested |
| `.d2` | `correlation_dimension_curve` | Implemented and tested |
| `bdsTest` / `bdstest_main` | `bds_test` | Implemented and tested |
| `wnnTest` | `white_neural_test` | Implemented and tested |
| `tnnTest` | `terasvirta_neural_test` | Implemented and tested |
| `runsTest` | `runs_test` | Implemented and tested |
| `tsTest` | `ts_test` | Implemented and tested |
| Plot methods and `doplot` branches | none | Excluded as plotting |
| S4 `fHTEST`, `ts`, and `timeSeries` wrappers | plain derived types/arrays | R infrastructure excluded |

The original C neighbor kernels are represented by exact full-dimensional
searches with direct, dynamic box-index, and automatic modes. The BDS pair
counts use a simpler direct algorithm while retaining the original statistic.
