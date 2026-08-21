# API coverage

| R function | Fortran counterpart | Status / notes |
|---|---|---|
| `dbp` | `dbp`, `dbp_scalar` | Implemented. Log PMF is the default conceptual form; scalar and vector APIs provided. |
| `rbp` | `rbp` | Implemented via independent Poisson components. |
| `bp.mle2` | `bp_mle2` | Implemented with one-dimensional profile maximization. |
| `bp.mle` | `bp_mle` | Implemented, including correlation, LLR/Wald p-values, CIs, observed and asymptotic variance calculations. |
| `lambda3.profile` | `lambda3_profile` | Numerical profile/grid and CI implemented; plotting omitted. |
| `bp.gof` | `bp_gof` | Implemented parametric-bootstrap dispersion test. Runtime reporting is omitted. |
| `bp.gof2` | `bp_gof2` | Implemented as compiled-loop equivalent of `bp_gof`; no R-style vectorization is needed. |
| `bp.contour` | `bp_probability_grid` | Probability-grid computation retained; contour/scatter graphics omitted. |
| R `table(x1,x2)` output | `make_bp_table` | Typed Fortran contingency table supplied. |

R-specific list names, matrix dimnames, plotting calls, `proc.time()`, and namespace/S3 presentation behavior are intentionally not reproduced.
