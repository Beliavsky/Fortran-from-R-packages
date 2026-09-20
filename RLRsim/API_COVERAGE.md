# API coverage

Coverage basis: distinct exported computational R functions in upstream `NAMESPACE`.

| R function | Status | Fortran mapping | Main compatibility difference |
|---|---|---|---|
| `LRTSim` | substantial | `rlrsim_api:lrt_sim` | No R parallel backend or exact R RNG stream. |
| `RLRTSim` | substantial | `rlrsim_api:rlrt_sim` | No R parallel backend or exact R RNG stream; numerical approximation branches are translated. |
| `exactLRT` | substantial | `rlrsim_api:exact_lrt_from_design` | Accepts explicit matrices and observed statistic instead of R fitted objects. |
| `exactRLRT` | substantial | `rlrsim_api:exact_rlrt_from_design` | Accepts explicit matrices and observed statistic instead of R fitted objects. |
| `extract.lmeDesign` | untranslated | none | Depends on R `nlme` object internals, formulas, factors, and model matrices. |

The unexported helper `extract.lmerModDesign` is not counted under the requested coverage convention.
