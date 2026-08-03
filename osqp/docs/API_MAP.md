# API map

| R package interface | Fortran interface | Notes |
|---|---|---|
| `solve_osqp(P,q,A,l,u,pars)` | `solve_osqp(q,solution,...,p,a,l,u,settings)` | One-shot solve. |
| `osqp(P,q,A,l,u,pars)` | `osqp_model_from_dense` or `osqp_model_from_sparse`, then `osqp_setup` | Creates a persistent solver. |
| `model@Solve()` | `osqp_solve_solver` | Returns `osqp_solution`. |
| `model@Update(q,l,u,Px,Px_idx,Ax,Ax_idx)` | `osqp_update` | Fortran indices are one-based. Sparsity patterns are preserved. |
| `model@WarmStart(x,y)` | `osqp_warm_start` | Either vector may be omitted. |
| `model@ColdStart()` | `osqp_cold_start` | Resets the ADMM iterate. |
| `model@GetParams()` | `osqp_get_settings` | Returns the active settings. |
| `model@GetDims()` | `osqp_get_dimensions` | Returns `n` and `m`. |
| `model@UpdateSettings()` | `osqp_update_settings` | OSQP setup-only settings remain unchanged, matching the C API. |
| `model@GetData()` | `solver%model` | The Fortran solver retains the unscaled user-level model. |
| `osqpSettings()` | `osqp_settings()` | Defaults match the R package/OSQP 1.0 defaults. |
| Dense/`dgCMatrix`/triplet coercion | `osqp_csc_from_dense`, `osqp_csc_from_triplet` | Duplicate triplets are summed. |

## Result mapping

`osqp_solution` contains `x`, `y`, `prim_inf_cert`, `dual_inf_cert`, the text status, and an `osqp_info` record. The information record includes objective values, residuals, duality gap, iterations, rho updates/estimate, timing, primal-dual integral, and relative KKT error.

## Omitted R infrastructure

S7/S3 objects, deprecated `$` method dispatch, `Matrix` and `slam` classes, R vector recycling, names, printing methods, CLI-formatted errors, vignettes, and R test infrastructure are not translated.
