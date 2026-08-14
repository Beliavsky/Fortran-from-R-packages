# API mapping

| R / C++ routine | Fortran counterpart | Notes |
|---|---|---|
| `mixsqp` | `fit_mixsqp` | Full preprocessing, EM + SQP and postprocessing. |
| `mixsqp_control_default` | `mixsqp_default_control` / `mixsqp_control` | Same numerical defaults; `maxiter_activeset=0` means automatic `min(20,m+1)`. |
| `mixobjective` | `mixobjective` | Normalizes weights and mixture vector before evaluation. |
| `mixobj` / `compute_objective` | internal `mixsqp_utils::mixobjective` | Inner objective used without renormalizing SQP iterates. |
| `mixem_rcpp` | `run_mixem` | Multi-step helper used internally. |
| `mixem_update` | `mixem_update` | Same stabilized posterior normalization. |
| `compute_grad` | `compute_grad_hessian` | Full or factorized likelihood path. |
| `activesetqp` | `active_set_qp` | Nocedal-Wright active-set method. |
| `compute_activeset_searchdir` | `compute_searchdir` | Cholesky with added multiple of identity. |
| `backtracking_line_search` | `line_search` | Same sufficient-decrease test and feasibility handling. |
| `tsvd` | `truncated_svd` | LAPACK `DGESDD` backend; same singular-value threshold semantics. |
| `simulatemixdata` | `simulate_mix_data` | Normal / normal+t examples and likelihood matrix construction. |
| `normalize.likelihoods` | `normalize_likelihoods` | Row maximum normalization. |
| `normalize.loglikelihoods` | `normalize_loglikelihoods` | Exponentiate after rowwise max subtraction. |
| `normalize.rows` | `normalize_rows_with_logscale` | Internal row normalization with log scaling vector. |
| `logspace` | `logspace` | Internal helper. |

R input validation is represented by Fortran shape/value checks and `error
stop` for invalid inputs. R debugging options, console formatting, names,
data frames and serialization are interface/UI facilities and are not ported.
