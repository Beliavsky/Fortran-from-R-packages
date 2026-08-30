# Computational coverage

| Upstream facility | Native Fortran coverage |
|---|---|
| `parm` | `make_parm` |
| `modelparm` numerical result | `parm_type`, `make_parm`, `parm_from_iid` |
| `glht` matrix hypotheses | `glht_fit` |
| `glht` default identity | `glht_identity` |
| `cftest` | `glht_coefficients` + `glht_test` |
| `univariate` | `glht_test(...,'univariate',...)` |
| `adjusted('single-step')` | native max-t normal/t integration |
| `adjusted('free')` | native free step-down |
| `adjusted('Shaffer')` | native closed testing + `maxsets` |
| `adjusted('Westfall')` | native closed testing + correlated local max tests |
| R `p.adjust` families | none, Bonferroni, Holm, Hochberg, Hommel, BH/FDR, BY |
| `Ftest`, `Chisqtest` | `glht_global_test` |
| adjusted/univariate `calpha` | `glht_critical_value` |
| `confint.glht` | `glht_confint` |
| `contrMat` | all ten documented contrast families |
| `cld` numerical algorithm | `compact_letter_display` and Tukey convenience wrappers |
| `mmm` covariance | `mmm_parm_from_iid` |
| `mlf` block assembly | `block_diagonal_matrix` |
| `mcp` contrast families | `contr_mat`; caller supplies numeric design mapping |
| character/expression hypothesis parser | omitted as R-language parsing |
| model-specific coefficient extraction | omitted as R object-system integration |
| print/plot/S3 methods | omitted as presentation infrastructure |
