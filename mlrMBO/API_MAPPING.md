# API mapping

This table maps important non-plotting mlrMBO computations to the Fortran
port. R S3/list constructors are represented by derived types and ordinary
procedure calls rather than duplicated literally.

| Upstream mlrMBO | Fortran counterpart | Status |
|---|---|---|
| `mbo`, `initSMBO`, `updateSMBO`, `finalizeSMBO` | `mbo`, `mbo_continue`, `finalize_mbo`, `mbo_result` | Implemented computational workflow |
| `makeMBOControl`, `setMBOControl*` | `mbo_control`, `init_control`, public fields | Implemented |
| `OptState`, `OptProblem`, `OptResult`, opt path | `mbo_space`, `mbo_path`, `mbo_result` | Implemented numerically |
| `crit.mr` / mean response | `crit_mean`, `eval_single_criterion` | Implemented |
| `crit.se` | `crit_se` | Implemented |
| `crit.ei` | `crit_ei`, `expected_improvement` | Implemented |
| `crit.cb`, `crit.cb1`, `crit.cb2` | `crit_cb`, `cb_lambda` | Implemented |
| adaptive CB | `crit_adacb`, `cb_lambda_start/end` | Implemented |
| `crit.aei` | `crit_aei` | Implemented |
| `crit.eqi` | `crit_eqi` | Implemented |
| DIB SMS/epsilon C kernels | `sms_indicator_values`, `eps_indicator_values` | Implemented natively |
| `getNonDominatedPoints`, `isDominated` | `nondominated_points`, `dominated_mask` | Implemented |
| dominated hypervolume / contributions | `dominated_hypervolume`, `hypervolume_contributions` | Implemented exact recursive algorithm |
| `getMultiObjRefPoint` | `reference_point` | Implemented (`all`, `front`, `const`) |
| ParEGO | `mo_parego`, `parego_weights`, `parego_scalarize` | Implemented |
| DIB | `mo_dib` | Implemented, including sequential batch front augmentation |
| MSPOT | `mo_mspot` | Implemented with native candidate-pool Pareto search instead of external NSGA-II |
| constant liar | `batch_cl` | Implemented |
| parallel confidence bounds | `batch_cb` | Implemented sequentially; lambda distribution preserved |
| random interleaving | `interleave_random_points` | Implemented |
| focus search | `focus_search` | Implemented |
| proposed-point filtering | `filter_proposed`, `filter_tol` | Implemented |
| `trafoLog`, `trafoSqrt` | `trafo_log`, `trafo_sqrt` and inverses | Implemented |
| target/evaluation budgets | `max_evals`, `max_iter`, `target_value` | Implemented |
| `mboContinue` | `mbo_continue` | Implemented |
| `makeMBOLearner`, arbitrary mlr learners | DiceKriging `km_model` surrogate | Framework substitution; external learner adapters deferred |
| CMA-ES/rgenoud/EA infill adapters | - | External-package adapters deferred |
| MOI-MBO | reserved `batch_moimbo` | Deferred; explicit error if selected |
| plotting/example rendering | - | Intentionally omitted |
| R parallel/data.table/ParamHelpers plumbing | - | Intentionally omitted |

## Parameter representation

`mbo_space` stores one encoded real column per parameter. Integer parameters
are rounded; categorical parameters use integer level codes stored as reals.
A simple conditional parameter may specify a categorical parent and the
required parent level. Inactive values are canonicalized to the parameter's
lower bound. The objective callback receives these repaired encoded values.
