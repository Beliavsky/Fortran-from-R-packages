# API mapping

The table maps the exported DiceDesign 1.10 numerical API to the Fortran API.
R-only plotting/display arguments are omitted.

| DiceDesign R export | Fortran counterpart | Notes |
|---|---|---|
| `coverage` | `coverage` | coefficient of variation of nearest-neighbor distances |
| `discrepancyCriteria` | `discrepancy_all`, `discrepancy_value` | seven discrepancy families |
| `discrepESE_LHS` | `discrep_ese_lhs` | ESE; full objective recomputation |
| `discrepSA_LHS` | `discrep_sa_lhs` | GEOM, LINEAR, GEOM_MORRIS and MC profiles |
| `dmaxDesign` | `dmax_design` | D-max stochastic design |
| `factDesign` | `fact_design` | full factorial design |
| `faureprimeDesign` | `faureprime_design` | prime-dimension Faure construction |
| `lhsDesign` | `lhs_design` | randomized or centered LHS |
| `maximinESE_LHS` | `maximin_ese_lhs` | minimizes phi-p |
| `maximinSA_LHS` | `maximin_sa_lhs` | SA maximin/phi-p optimizer |
| `meshRatio` | `mesh_ratio` | nearest-neighbor mesh ratio |
| `mindist` | `mindist` | minimum pairwise Euclidean distance |
| `mstCriteria` | `mst_criteria` | Prim MST + edge mean/sd; plotting omitted |
| `nolhDesign` | `nolh_design` | embedded upstream NOLH tables |
| `nolhdrDesign` | `nolhdr_design` | embedded upstream NOLHDR tables |
| `olhDesign` | `olh_design` | recursive orthogonal LHS |
| `phiP` | `phi_p` | phi-p distance criterion |
| `rss2d` | `rss2d` | numerical RSS result; graphics omitted |
| `rss3d` | `rss3d` | numerical RSS result; graphics omitted |
| `runif.faure` | `runif_faure` | Faure low-discrepancy sequence |
| `scaleDesign` | `scale_design` | linear or empirical-CDF scaling |
| `straussDesign` | `strauss_design` | explicit local RNG state |
| `unif.test.quantile` | `unif_test_quantile` | supported upstream quantiles |
| `unif.test.statistic` | `unif_test_statistic` | optional spacing transform |
| `unscaleDesign` | `unscale_design` | linear or empirical-quantile inverse |
| `wspDesign` | `wsp_design` | WSP thinning with selected indices |
| `xDRDN` | `xdrdn_transform` | numeric rescale/round transform only |

Internal R helper functions used solely to maintain incremental LHS objective
updates are not exposed one-for-one. The Fortran SA/ESE routines evaluate the
same public objective directly after each elementary proposal, which is simpler
and avoids unstable divisions in boundary cases.
