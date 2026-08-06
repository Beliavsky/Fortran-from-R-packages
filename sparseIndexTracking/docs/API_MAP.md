# API map

| Upstream R symbol | Modern Fortran symbol | Notes |
|---|---|---|
| `spIndexTrack` | `fit_sparse_index_tracking` | Typed result API. |
| `spIndexTrack` | `sp_index_track` / `spIndexTrack` | Direct weight-vector compatibility API. |
| `eteMMupdate` | internal `mm_update`, `measure_ete` | ETE MM update. |
| `drMMupdate` | internal `mm_update`, `measure_dr` | Downside-risk MM update. |
| `heteMMupdate` | internal `mm_update`, `measure_hete` | Symmetric Huber MM update. |
| `hdrMMupdate` | internal `mm_update`, `measure_hdr` | One-sided Huber MM update. |
| `bisection` | `bisection` | Solves the capped-simplex KKT equation. |
| objective expressions | `tracking_objective` | Supports all four measures and optional source HDR behavior. |

## Public result fields

`type(sparse_index_fit)` contains:

- `weights`
- `objective`
- `iterations`
- `outer_iterations`
- `cardinality`
- `info`
- `converged`
- `message`

## Status codes

- `sit_success`
- `sit_invalid_argument`
- `sit_dimension_error`
- `sit_infeasible_bounds`
- `sit_degenerate_data`
- `sit_iteration_limit`
- `sit_numerical_error`
