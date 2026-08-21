# API coverage

| R API | Fortran API | Status / notes |
|---|---|---|
| `paras()` | `paras`, `paras_from_values`, `paras_from_p_theta` | Implemented with `paras_type`. |
| `p()` / `p<-` | `p`, `set_p` | Implemented. `set_p` accepts either `k-1` or `k` entries. |
| `theta()` / `theta<-` | `theta`, `set_theta` | Implemented; scalar and matrix replacement supported. |
| `pnames()` | `paras_type%pnames`, `mb_type%pnames` | Stored as fixed-length strings; no S4 accessor required. |
| `getVals()` | `paras_type%vals` | Direct typed access. |
| `length(paras)` | `paras_dimension` | Implemented. |
| `MB()` | `make_mb` | Implemented as `mb_type`; validates `0 <= counts <= m`. |
| `counts()` / `getM()` | `mb_type%counts`, `mb_type%m` | Direct typed access. |
| `lmultinomial()` | `lmultinomial` | Implemented. |
| `multinomial()` | `multinomial` | Implemented. |
| `MM_single()` | `mm_single`, `mm_single_log` | Implemented. |
| `dMM()` | `dmm` | Implemented. |
| `NormC()` | `normc`, `normc_log` | Implemented with exact composition enumeration + log-sum-exp. |
| `MM()` | `mm_loglik` | Implemented dispatcher. |
| `MM_allsamesum()` | `mm_allsamesum` | Implemented. |
| `MM_differsums()` | `mm_differsums` | Implemented. |
| `MM_allsamesum_A()` | `mm_allsamesum_a` | Implemented with upstream formula, including historical `NormC` term. |
| `MM_differsums_A()` | `mm_differsums_a` | Implemented. |
| `MM_support()` | `mm_support` | Implemented. |
| `optimizer()` | `optimizer` | Implemented. Default BFGS analogue of `nlm`; optional Nelder-Mead. |
| `optimizer_allsamesum()` | `optimizer_allsamesum` | Implemented. |
| `optimizer_differsums()` | `optimizer_differsums` | Implemented. |
| `gunter(matrix/data.frame)` | `gunter` | Implemented as `gunter_type`. |
| `gunter(MB)` | `gunter_mb` | Implemented as `gunter_mb_type`. |
| `as.array.*`, Oarray conversion | -- | Omitted as R/Oarray container machinery; support table is available directly. |
| `Lindsey()` | `lindsey`, `lindsey_fit` | Implemented with native Poisson IRLS/Newton fit. |
| `Lindsey_MB()` | `lindsey_mb` | Implemented for the upstream bivariate case. Returns `glm_fit_type`. |
| `rMM()` | `rmm` | Implemented Metropolis-Hastings sampler. |
| `suffstats()` | `suffstats` | Implemented. |
| `summary.suffstats()` | divide fields by `nobs` | Representation/printing convenience not duplicated as a generic method. |
| `expected_suffstats()` | `expected_suffstats` | Implemented. |
| print/show methods | -- | Omitted. |

## Dependency mappings

- R `partitions::compositions` -> supplied `partitions-fortran::compositions`.
- R `quadform::quad.form` / `quad.tdiag` / `quad.ttrace` -> supplied
  `quadform-fortran` operations.  The port uses `quad_tdiag` for row-wise
  interaction terms; this is algebraically the same computation.
- R `stats::glm(..., family=poisson)` -> native `poisson_glm_fit`.
- R `stats::nlm` -> native finite-difference BFGS smooth optimizer.
- R `stats::optim(..., method=Nelder-Mead)` -> native Nelder-Mead optimizer.
