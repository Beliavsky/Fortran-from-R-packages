# API map

Upstream: GPArotation 2026.8-2.

## Gradient projection and random starts

| R | Fortran |
|---|---|
| `GPForth` | `gpa_rotation:gpforth` |
| `GPFoblq` | `gpa_rotation:gpfoblq` |
| `GPFRSorth` | `gpa_api:gpfrsorth` |
| `GPFRSoblq` | `gpa_api:gpfrsoblq` |
| `Random.Start` | `gpa_api:random_start` |

`rotation_options` exposes convergence tolerance, maximum iterations,
`legacy`/`bb`/`cayley` step choice (Cayley orthogonal only), non-monotone line
search window, normalization mode, and random-start count.

## Rotation wrappers

| R | Fortran |
|---|---|
| `oblimin` | `oblimin` |
| `quartimin` | `quartimin` |
| `targetT`, `targetQ` | `target_t`, `target_q` |
| `pstT`, `pstQ` | `pst_t`, `pst_q` |
| `oblimax` | `oblimax` |
| `binormamin` | `binormamin` |
| `entropy` | `entropy` |
| `quartimax` | `quartimax` |
| `Varimax` | `varimax` |
| `simplimax` | `simplimax` |
| `bentlerT`, `bentlerQ` | `bentler_t`, `bentler_q` |
| `tandemI`, `tandemII` | `tandem_i`, `tandem_ii` |
| `geominT`, `geominQ` | `geomin_t`, `geomin_q` |
| `bigeominT`, `bigeominQ` | `bigeomin_t`, `bigeomin_q` |
| `cfT`, `cfQ` | `cf_t`, `cf_q` |
| `infomaxT`, `infomaxQ` | `infomax_t`, `infomax_q` |
| `mccammon` | `mccammon` |
| `bifactorT`, `bifactorQ` | `bifactor_t`, `bifactor_q` |
| `equamax` | `equamax` |
| `parsimax` | `parsimax` |
| `varimin` | `varimin` |
| `lpT`, `GPForth.lp` | `lp_t`, `lp_rotate(...,.true.,...)` |
| `lpQ`, `GPFoblq.lp` | `lp_q`, `lp_rotate(...,.false.,...)` |
| `eiv` | `eiv` / `eiv_rotate` |
| `echelon` | `echelon` / `echelon_rotate` |

## Criterion functions

Every upstream smooth/non-smooth `vgQ.*` criterion is directly represented in
`gpa_criteria`:

`vgq_oblimin`, `vgq_quartimin`, `vgq_cf`, `vgq_target`, `vgq_pst`,
`vgq_entropy`, `vgq_infomax`, `vgq_mccammon`, `vgq_geomin`,
`vgq_simplimax`, `vgq_bifactor`, `vgq_bigeomin`, `vgq_tandem1`,
`vgq_tandem2`, `vgq_oblimax`, `vgq_bentler`, `vgq_quartimax`,
`vgq_varimax`, `vgq_binormamin`, `vgq_varimin`, and `vgq_lp_wls`.

`evaluate_criterion` dispatches by the upstream criterion name.

## Normalization and diagnostics

- `NormalizingWeight`: `normalizing_weights` (none, Kaiser, Cureton-Mulaik).
- `residuals.GPArotation`: `residual_matrix`.
- `calc_fitstats`: `calc_fitstats` (df, chi-square, SRMR, RMSEA and CI).
- `calc_FSI`: `calc_fsi`.
- `calc_AUC`: `calc_auc`.
- `calc_simplicity`: `calc_simplicity`.
- `calc_hyperplane`: `calc_hyperplane`.

## Not translated

Graphics and R presentation infrastructure are excluded:
`plot.GPArotation`, trajectory/landscape plots, `plot2fOrthComparison`, and S3
printing/summary formatting. `factanal` extraction is replaced by direct matrix
arguments plus explicit correlation/sample-size inputs to diagnostics.
