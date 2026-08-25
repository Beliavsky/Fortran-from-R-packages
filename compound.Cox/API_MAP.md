# API map

| Upstream R function | Fortran procedure | Status |
|---|---|---|
| `CG.Clayton` | `cg_clayton` | implemented |
| `CG.Frank` | `cg_frank` | implemented |
| `CG.Gumbel` | `cg_gumbel` | implemented |
| `CG.test` | `cg_test` | implemented, non-plotting |
| `X.pathway` | `x_pathway` | implemented |
| `X.tag` | `x_tag` | implemented |
| `cindex.CV` | `cindex_cv` | implemented |
| `compound.reg` | `compound_reg` | implemented, including variance correction |
| `dependCox.reg` | `depend_cox_reg` | implemented, including censor coefficient/baselines |
| `dependCox.reg.CV` | `depend_cox_reg_cv` | implemented |
| `surv.factorial` | `surv_factorial` | implemented, including jackknife/F test |
| `uni.Wald` | `uni_wald` | implemented |
| `uni.score` | `uni_score` | implemented |
| `uni.selection` | `uni_selection` | implemented, including CVL/FDR permutation |

## Presentation-only differences

R plotting options (`S.plot`, `CC.plot`, plotting in `dependCox.reg.CV`, and
`surv.factorial`) are not part of the Fortran API. R list objects are represented
by typed derived results. `CG.test` and `surv_factorial` select one of the three
built-in copulas by name instead of accepting an arbitrary R function object.
