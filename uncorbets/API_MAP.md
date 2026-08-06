# API map

| Upstream R routine | Fortran routine | Notes |
|---|---|---|
| `torsion(..., model="pca")` | `torsion`, `torsion_pca` | Symmetric Jacobi eigendecomposition; descending eigenvalues and upstream sign convention. |
| `torsion(..., model="minimum-torsion", method="approximate")` | `torsion`, `torsion_minimum` | Direct translation of the Riccati-root approximation. |
| `torsion(..., model="minimum-torsion", method="exact")` | `torsion`, `torsion_minimum` | Direct translation of the iterative minimum-torsion algorithm. |
| `effective_bets` | `effective_bets` | Returns `effective_bets_result` with probability vector and ENB. |
| `max_effective_bets` | `max_effective_bets` | Same simplex constraints; independent projected-gradient optimizer replaces `NlcOptim::solnl`. |
| `sqrtm` | `sqrtm` | Symmetric positive-semidefinite matrix square root. |
| `is_quadratic` | Input validation | Square-matrix checks are built into public procedures. |
| `is_col_named` | Omitted | Matrix names are R metadata and have no native Fortran equivalent. |
