# API coverage

Coverage basis: distinct exported computational R functions and registered S3 computational methods in upstream `NAMESPACE`. Presentation-only `print.poLCA` and `plot.poLCA` are excluded. Internal helpers and native `.C` routines are not separate denominator entries, although their numerical operations are translated and exposed where useful.

Package status: **substantial**. Coverage: **10 of 10 (100%)**.

The 100% mapping fraction does not mean complete R compatibility. The principal intentional differences are numeric-array APIs instead of formulas/data frames/S3 objects, no plotting or formatted printing, and no attempt to reproduce R's RNG stream bit-for-bit.

| R function | Fortran procedure(s) | Status | Main differences |
|---|---|---|---|
| `poLCA` | `polca_fit` | substantial | Formula/model-frame handling, R object construction, console restart messages, and plotting are omitted. Numerical LCA/LCR fitting, missing responses, starts, fit statistics, and standard errors are implemented. |
| `poLCA.simdata` | `polca_simdata` | substantial | Uses numeric arrays; missingness is independent cell-wise Bernoulli rather than R's sampled index-pair mechanism; RNG stream differs. |
| `poLCA.reorder` | `polca_reorder_probs` | complete | Rectangular rank-3 Fortran probability array replaces an R list of matrices. |
| `poLCA.table` | `polca_table_1d`, `polca_table_2d` | substantial | Item indices/categories replace formulas, names, factors, and presentation labels. |
| `poLCA.entropy` | `polca_entropy` | complete | Same full-cell probability entropy. |
| `poLCA.predcell` | `polca_predcell` | complete | Same complete-cell probability calculation with numeric validation. |
| `poLCA.posterior` | `polca_posterior` | substantial | Numeric matrices replace R coercion/model reconstruction. |
| `rmulti` | `polca_rmulti` | substantial | Uniform draws are explicit in the low-level routine for deterministic testing. |
| `coef.poLCA` | `polca_coef` | complete | Zero-size matrix represents absence of regression coefficients instead of R `NA`. |
| `vcov.poLCA` | `polca_vcov` | complete | Zero-size matrix represents absence of coefficient covariance instead of R `NA`. |

## Internal upstream numerical operations translated

The translation also contains direct Fortran equivalents of the computational work in `poLCA.ylik.C`, `poLCA.postClass.C`, `poLCA.probHat.C`, `poLCA.dLL2dBeta.C`, `poLCA.updatePrior`, `poLCA.se`, and complete-cell compression/goodness-of-fit calculations. These are internal upstream helpers and therefore are not counted as exported-function coverage entries.
