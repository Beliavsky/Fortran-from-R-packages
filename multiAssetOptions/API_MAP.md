# API map

| Upstream R function | Fortran equivalent | Notes |
|---|---|---|
| `nodeSpacer` | `node_spacer` | Uniform/nonuniform spacing and mesh shifts 0, 1, and 2. |
| `payoff` | `payoff_values` | Returns a classically ordered rank-one state vector rather than an R array. |
| `matrixFDM` | `build_fdm_operator` | Returns native `csr_matrix`; uses equivalent tensor-product derivative stencils. |
| `multiAssetOption` | `price_multi_asset` | Uses typed `pricing_config` and returns `pricing_result`. |
| `plotOptionValues` | Omitted | Plotting/animation is outside the computational scope. |

Additional Fortran helpers include `initialize_config`, `build_grid`,
`interpolate_value`, `csr_matvec`, `csr_to_dense`, and index conversion
routines.
