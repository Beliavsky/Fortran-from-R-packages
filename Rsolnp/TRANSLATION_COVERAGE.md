# Translation coverage

| R export | Fortran entry point | Status | Notes |
|---|---|---|---|
| `solnp` | `solnp` | Complete | Uses the common translated solver |
| `csolnp` | `csolnp` | Complete | Analytic derivatives optional |
| `csolnp_ms` | `csolnp_ms` | Complete | Sequential deterministic multistart |
| `gosolnp` | `gosolnp` | Complete | Alias of translated multistart driver |
| `startpars` | `startpars` | Complete | Halton starts plus feasibility improvement |
| `kkt_diagnose` | `kkt_diagnose` | Complete | Explicit diagnostics derived type |
| `solnp_standardize_problem` | same | Complete | Typed standard-form metadata |
| `solnp_problem_suite` | same | Partial benchmark definitions | 18 executable definitions; full source retained |
| `solnp_problems_table` | same | Complete | 77 registry entries, including HS65 |

## Computational internals

Translated/adapted internals include augmented-objective evaluation, projected
BFGS, Armijo search, feasibility restoration, finite-difference gradients and
Jacobians, bounded slack handling, penalty/multiplier updates, regularized KKT
multiplier estimation, start ranking, and benchmark callbacks.

## Out of scope

- R S3 print and summary methods
- R list/formula/ellipsis dispatch
- Parallel execution through `parallel` or `future.apply`
- Vignettes and R package build infrastructure
- Exact R RNG replication for benchmark data
- Executable ports of every benchmark row not needed by the solver API
