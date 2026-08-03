# Translation coverage

ICSNP 1.1-3 exports 20 computational entry points. All 20 are represented.

| R export | Fortran procedure | Status |
|---|---|---|
| `rank.ctest` | `rank_ctest`, `rank_ctest_groups` | translated |
| `rank.ictest` | `rank_ictest` | translated |
| `HotellingsT2` | `HotellingsT2` | translated |
| `spatial.median` | `spatial_median` | translated |
| `spatial.sign` | `spatial_sign` | translated |
| `HP1.shape` | `HP1_shape` | translated |
| `HR.Mest` | `HR_Mest` | translated |
| `duembgen.shape` | `duembgen_shape` | translated |
| `tyler.shape` | `tyler_shape` | translated |
| `ind.ictest` | `ind_ictest` | translated; self-contained FOBI transform |
| `ind.ctest` | `ind_ctest` | translated |
| `pair.prod` | `pair_prod` | translated |
| `pair.sum` | `pair_sum` | translated |
| `pair.diff` | `pair_diff` | translated |
| `duembgen.shape.wt` | `duembgen_shape_wt` | translated |
| `symm.huber` | `symm_huber` | translated |
| `symm.huber.wt` | `symm_huber_wt` | translated |
| `hl.loc` | `hl_loc` | translated |
| `vdw.loc` | `vdw_loc` | translated |
| `HP.loc.test` | `HP_loc_test` | translated |

Internal C++ spatial-rank and signed-rank kernels are also exposed as
`spatial_ranks` and `signed_ranks`.

Formula/S3/ICS-object wrappers and R package infrastructure are documented as
out of scope in `original/OMITTED.md`.
