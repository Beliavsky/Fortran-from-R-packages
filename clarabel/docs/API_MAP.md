# API map

## R to Fortran

| R interface | Fortran counterpart | Notes |
|---|---|---|
| `clarabel()` | `clarabel_solve_problem()` | One-shot solve. |
| `clarabel_solver()` / `ClarabelSolver$new()` | `clarabel_solver_type%initialize()` | Persistent solver object. |
| `solver_solve()` | `clarabel_solver_type%solve()` | Returns `clarabel_solution`. |
| `solver_update()` | `clarabel_solver_type%update()` | Optional `P`, `A`, `q`, and `b`; matrix patterns must match. |
| `solver_is_update_allowed()` | `clarabel_solver_type%is_update_allowed()` | Reflects presolve, zero dropping, and chordal decomposition. |
| `clarabel_control()` | `clarabel_settings` / `default_clarabel_settings()` | Typed fields replace an R list. |
| `make_csc_matrix()` | `csc_from_dense()`, `csc_from_arrays()`, `csc_from_triplets()` | Zero-based CSC is stored internally. |
| `make_csc_symm_matrix()` | `csc_from_symmetric_upper()` | Stores only the upper triangle of `P`. |
| solver status descriptions | `status_name()` and status constants | Includes callback termination from the Rust core. |

## Cone specification

| R cone entry | Fortran constructor |
|---|---|
| `z` | `zero_cone(dim)` |
| `l` | `nonnegative_cone(dim)` |
| `q` | `second_order_cone(dim)` |
| `ep` | `exponential_cone()` once per exponential cone |
| `p` | `power_cone(alpha)` |
| `gp = list(a=..., n=...)` | `generalized_power_cone(alpha, dim)` |
| `s` | `psd_triangle_cone(matrix_order)` |

The total cone dimension must equal `size(A,1)`. A PSD cone of matrix order `n` contributes `n*(n+1)/2` rows. A generalized-power cone contributes `size(alpha)+dim` rows.

## Settings codes

`direct_solve_method` uses:

- `direct_solver_auto`
- `direct_solver_qdldl`
- `direct_solver_faer`
- `direct_solver_mkl`
- `direct_solver_panua`

Availability depends on the features used to build Clarabel.rs. The bundled default bridge enables the standard solver plus SDP support; QDLDL is the portable default.

`chordal_decomposition_merge_method` uses:

- `chordal_merge_none`
- `chordal_merge_parent_child`
- `chordal_merge_clique_graph`

## Result mapping

`clarabel_solution` contains `x`, `z`, `s`, primal and dual objectives, solve time, iteration count, residuals, status, and a nested `clarabel_info` record. The latter exposes path parameters, costs, residuals, gaps, kappa/tau ratio, and linear-solver statistics.
