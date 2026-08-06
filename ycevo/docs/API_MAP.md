# R-to-Fortran API map

| R / Rcpp routine | Fortran routine or type | Notes |
|---|---|---|
| `epaker` | `epanechnikov` | Elemental function. |
| `qepa` | `epanechnikov_quantile` | Elemental inverse CDF. |
| `repa` | `epanechnikov_random` | Uses Fortran intrinsic RNG. |
| `calc_epaker_weights` | `calc_epaker_weights` | Allocatable matrix result. |
| `get_weights` | `get_weights` | Supports arbitrary units. |
| `range_idx_nonzero` | `range_nonzero` | Returns first/last indices separately. |
| `span2h` | `span_to_bandwidth` | Discrete bisection matching the support count. |
| `num_points_mat` | `count_maturing_bonds` | Operates on `bond_panel_t`. |
| `create_tau_ht` | `adjust_tau_grid` | Trims the sparse tail, removes interior sparse points, and widens gap-edge bandwidths. |
| `nelson_siegel` | `nelson_siegel` | Elemental function. |
| `get_yield_at`, `get_yield_at_vec` | `get_yield_at` | Elemental and therefore naturally vectorizable. |
| `generate_yield` | `generate_yield` | Allocates a maturity-by-date matrix. |
| `discount2yield` | `discount_to_yield` | Elemental function. |
| `get_cfp_slist` | `bond_panel_t` | Sparse lists replaced by a sorted flattened panel. |
| `calc_dbar_c`, `calc_dbar` | `calc_dbar` | Direct loops over nonzero cash-flow rows. |
| `calc_hhat_num_c`, `calc_hhat_num` | `calc_hhat_numerator` | Grouped by quotation day and bond identifier. |
| interpolation matrix in `estimate_yield` | `interpolation_weights` | Preserves the R package's weight orientation. |
| `estimate_yield` | `estimate_yield` | Single date/covariate point. |
| `ycevo` | `estimate_yield_surface` plus grid utilities | Explicit typed API instead of tibble/S3 orchestration. |
| `augment.ycevo`, `predict.ycevo` | `loess_predict`, `predict_discount`, `predict_yield` | Local quadratic smoothing followed by clamped interpolation. |
| `interp1`, `interp2` | `linear_interpolate`, `bilinear_interpolate` | Scalar interpolation functions. |
| `ycevo_data` | `simulate_bond_panel` | Calendar-free computational analogue using integer day indices. |
| Data import/export | `read_bond_panel_csv`, `write_yield_curve_csv` | Native formatted I/O. |
| `autoplot`, `plot`, `vis_kernel` | Not translated | Plotting intentionally omitted. |
