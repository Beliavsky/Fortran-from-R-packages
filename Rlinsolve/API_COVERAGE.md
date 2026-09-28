# API coverage

Coverage basis: distinct exported computational R functions in `NAMESPACE`.
Presentation, package-attachment, and internal helper functions are excluded.

**Status: substantial — 12 of 12 (100%) mapped.**

The percentage measures mapped computational functions, not complete R object,
RNG, sparse-backend, warning/message, or full R object compatibility. Dense matrix-valued right-hand sides are supported by rank-2 overloads for every solver.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `aux.fisch` | `rlinsolve` | `aux_fisch` | complete |
| `lsolve.bicg` | `rlinsolve` | `lsolve_bicg` | substantial |
| `lsolve.bicgstab` | `rlinsolve` | `lsolve_bicgstab` | substantial |
| `lsolve.cg` | `rlinsolve` | `lsolve_cg` | substantial |
| `lsolve.cgs` | `rlinsolve` | `lsolve_cgs` | substantial |
| `lsolve.cheby` | `rlinsolve` | `lsolve_cheby` | substantial |
| `lsolve.gmres` | `rlinsolve` | `lsolve_gmres` | substantial |
| `lsolve.gs` | `rlinsolve` | `lsolve_gs` | substantial |
| `lsolve.jacobi` | `rlinsolve` | `lsolve_jacobi` | substantial |
| `lsolve.qmr` | `rlinsolve` | `lsolve_qmr` | substantial |
| `lsolve.sor` | `rlinsolve` | `lsolve_sor` | substantial |
| `lsolve.ssor` | `rlinsolve` | `lsolve_ssor` | substantial |
