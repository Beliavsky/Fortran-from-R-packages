# API coverage

Status: substantial. 26 of 26 (100%) exported computational functions are mapped;
this is not a claim of full R interface compatibility. No computational exports
remain untranslated. The authoritative machine-readable inventory is in
`fpm.toml`. All procedures below are public through `ineq_api`.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `ineq` | `ineq_indices` | `ineq` | substantial |
| `Gini` | `ineq_indices` | `gini` | substantial |
| `RS` | `ineq_indices` | `rs` | substantial |
| `Atkinson` | `ineq_indices` | `atkinson` | substantial |
| `Theil` | `ineq_indices` | `theil` | substantial |
| `Kolm` | `ineq_indices` | `kolm` | substantial |
| `var.coeff` | `ineq_indices` | `var_coeff` | substantial |
| `entropy` | `ineq_indices` | `entropy` | substantial |
| `conc` | `ineq_indices` | `conc` | substantial |
| `Herfindahl` | `ineq_indices` | `herfindahl` | substantial |
| `Rosenbluth` | `ineq_indices` | `rosenbluth` | substantial |
| `pov` | `ineq_indices` | `pov` | substantial |
| `Watts` | `ineq_indices` | `watts` | substantial |
| `Sen` | `ineq_indices` | `sen` | substantial |
| `SST` | `ineq_indices` | `sst` | substantial |
| `Foster` | `ineq_indices` | `foster` | substantial |
| `Lc` | `ineq_lorenz` | `lc` | substantial |
| `Lc.mehran` | `ineq_lorenz` | `lc_mehran` | substantial |
| `Lasym` | `ineq_lorenz` | `lasym` | substantial |
| `theorLc` | `ineq_lorenz` | `theor_lc` | substantial |
| `Lc.dagum` | `ineq_lorenz` | `lc_dagum` | substantial |
| `Lc.exp` | `ineq_lorenz` | `lc_exp` | substantial |
| `Lc.lognorm` | `ineq_lorenz` | `lc_lognorm` | substantial |
| `Lc.pareto` | `ineq_lorenz` | `lc_pareto` | substantial |
| `Lc.singh` | `ineq_lorenz` | `lc_singh` | substantial |
| `major` | `ineq_lorenz` | `major` | substantial |

## Exclusions

`Pen` is presentation-only (plotting). Registered `plot.Lc`, `lines.Lc`,
`plot.theorLc`, and `lines.theorLc`, plus dataset `Ilocos`, are not computational
exports under this coverage basis. R S3 object construction, graphics, warnings,
partial string matching, implicit coercion and argument recycling are omitted.
`theorLc` is included because its selected distribution is numerically evaluated.

See API.md for the exact domains and documented special-value differences.
