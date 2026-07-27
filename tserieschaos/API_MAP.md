# API map

## Exported R routines

| Original routine | Modern Fortran procedure | Status |
|---|---|---|
| `embedd` | `delay_embed`, `delay_embed_lags`, `delay_embed_matrix` | Implemented and tested |
| `C2` | `correlation_integral` | Implemented and tested |
| `d2` | `correlation_dimension_curve` | Implemented and tested |
| `mutual` | `average_mutual_information` | Implemented and tested |
| `false.nearest` | `false_nearest_fraction`, `false_nearest_curve` | Implemented with direct, box, and auto search; tested |
| `lyap_k` | `lyapunov_stretching` | Implemented with direct, box, and auto search; tested |
| `lyap` | `lyapunov_linear_fit` | Implemented and tested |
| `recurr` | `recurrence_distance_matrix` | Numerical matrix implemented and tested; plotting excluded |
| `stplot` | `space_time_separation` | Numerical isolines implemented and tested; plotting excluded |
| `sim.cont` | `integrate_rk4`, `simulate_observed` | Tested numerical analogue using RK4 instead of LSODA |
| `lorenz.syst` | `lorenz_rhs`, `simulate_lorenz` | Implemented and tested |
| `rossler.syst` | `rossler_rhs`, `simulate_rossler` | Implemented and tested |
| `duffing.syst` | `duffing_rhs`, `simulate_duffing` | Implemented and tested |

## Internal C routines

| Original C routine | Modern Fortran procedure | Status |
|---|---|---|
| `C2` | `correlation_integral` | Implemented |
| `d2` | `correlation_dimension_curve` | Implemented with direct threshold counting |
| `mutual` | `average_mutual_information` | Implemented, including the original one-marginal finite-sample formula |
| `falseNearest` | `false_nearest_fraction` | Implemented with selectable direct or box search |
| `find_knearests` | `find_k_nearests` | Implemented with selectable direct or box search |
| `follow_points` | `follow_neighbor_points` | Implemented |
| `stplot` | `space_time_separation` | Implemented with the original 1000-bin and 10-isoline construction |
| `init_boxSearch`, `fill_boxes` | Internal dynamic box-index construction | Implemented and tested |
| `find_nearests` | Internal exact radius query over neighboring boxes | Implemented and tested against direct search |

## Neighbor-search API

The following procedures accept `search_method="direct"`, `"box"`, or `"auto"`:

- `false_nearest_fraction`
- `false_nearest_curve`
- `find_k_nearests`
- `lyapunov_stretching`

The default is `auto`. The scalar and k-nearest routines also expose optional `distance_evaluations`; routines that make one search selection expose optional `method_used`.

The returned neighbors are ordered by ascending Euclidean distance and then ascending sample index. This makes tie handling deterministic and identical between direct and box search.

## Excluded infrastructure

The following are not computational omissions:

- `plot.ami`, `plot.d2`, `plot.false.nearest`
- `print.d2`, `print.false.nearest`
- `filled.contour`, `scatterplot3d`, and other graphics calls
- S3 classes and methods
- R `ts` metadata, date windows, frequencies, and `zoo`-style indexing
- R console progress messages
- packaged `.RData` examples

## Array and index conventions

- Fortran indices are one-based.
- Missing neighbor slots are returned as `-1`.
- `average_mutual_information` is allocated with lower bound zero, so `ami(0)` is lag zero.
- Delay embedding uses forward coordinates `x(i), x(i+d), ...`, matching the coordinates used by the original C algorithms.
- `lyapunov_linear_fit` accepts sample indices and an optional `dt`; the fitted slope is per sample when `dt` is omitted and per time unit when it is supplied.
