# API map

| R API | Fortran API | Coverage |
|---|---|---|
| `pcor(x, method)` | `pcor(x, result, method)` or `partial_correlation` | Direct numerical counterpart |
| `spcor(x, method)` | `spcor(x, result, method)` or `semi_partial_correlation` | Direct numerical counterpart |
| `pcor.test(x, y, z, method)` | generic `pcor_test(x, y, z, result, method)` | Vector or matrix `z` |
| `spcor.test(x, y, z, method)` | generic `spcor_test(x, y, z, result, method)` | Vector or matrix `z` |
| `method="pearson"` | `ppcor_pearson` | Complete |
| `method="spearman"` | `ppcor_spearman` | Complete, average tied ranks |
| `method="kendall"` | `ppcor_kendall` | Complete, Kendall tau-b |
| `MASS::ginv` fallback | native symmetric Moore-Penrose inverse | Complete numerical counterpart |
| R list/data-frame returns | `ppcor_result`, `ppcor_test_result` | Typed adaptation |

`method_from_name()` converts a character method name to the integer constants,
and `method_name()` performs the reverse conversion.
