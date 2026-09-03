# API coverage

This document distinguishes translated numerical algorithms from R-specific
adapters. Names in the first column refer to the supplied pbkrtest 0.5.5 source.

| Upstream area / routine | Fortran API | Coverage |
|---|---|---|
| `KR_Sigma_G2.R`: `.get_SigmaG` numerical Sigma/G construction | `build_sigma_g` | Core decomposition translated; caller supplies expanded random-effect `Z` blocks and covariance matrices instead of an `lmerMod` object. |
| `KR_vcovAdj.R`: `vcovAdj_internal` | `vcov_adjust_kr` | Numerical algorithm translated, including `P`, `Q`, information matrix, inverse/pseudoinverse branch, `W`, and adjusted covariance. |
| `KR_modcomp.R`: `.KR_adjust` | `kr_adjust` | F statistic, moment corrections, scaling, denominator df, and p-values translated. The upstream diagnostic convention for the unscaled statistic is preserved. |
| `get_info_functions.R`: `Lb_ddf` | `lb_ddf` | Translated. |
| `get_info_functions.R`: `ddf_Lb` | `ddf_lb_scalar` | Scalar-contrast numerical algorithm translated. |
| `SAT_modcomp.R`: `compute_auxiliary` | `compute_auxiliary_numeric` | Hessian eigen/pseudoinverse and covariance-Jacobian algorithm translated using callbacks and the shared `numDeriv`; `lmerMod` environment extraction is omitted. |
| `SAT_modcomp.R`: SAT F computation | `satterthwaite_test` | Translated from array inputs. |
| `SAT_modcomp.R`: `get_Fstat_ddf` | `get_fstat_ddf` | Translated, including equal-df and `nu <= 2` branches. |
| `PB_modcomp.R`: `PBcompute_p_values` | `bootstrap_p_values` | Calibration translated: direct PB, chi-square, Bartlett, gamma moment match, and moment-matched F. Reference simulation/refitting is omitted. |
| `get_info_functions.R`: model-specific `getLRT.*` numeric core | `likelihood_ratio_test` | Generic log-likelihood/parameter-count calculation translated. Model-object extraction is omitted. |
| `linear_algebra.R`: `compare_column_space` | `compare_column_space` | Translated. |
| `KR_utils.R`: `orthogonal_complement` | `orthogonal_complement` | Translated with full SVD instead of R QR; returns an equivalent orthonormal complement. |
| `model_coerce.R`: `force_full_rank` | `force_full_rank` | Translated with SVD; preserves row space rather than R's specific QR basis signs. |
| `model_coerce.R`: `make_model_matrix` | `make_model_matrix` | Translated. |
| `model_coerce.R`: `make_restriction_matrix` | `make_restriction_matrix` | Translated with SVD/projector linear algebra; the returned full-row-rank restriction matrix spans the same row space, though basis signs/rotations may differ from R QR output. |
| `PBrefdist`, `PBmodcomp`, `KRmodcomp`, `SATmodcomp` model refits | — | R/lme4 object orchestration omitted. Numerical kernels needed after model quantities are supplied are translated above. |
| Formula parsing, S3 dispatch/coercion, `broom`/`dplyr` output, printing/summary methods | — | R-specific interface code intentionally omitted. |
| Parallel cluster setup and R RNG/reference simulation | — | R-specific orchestration/RNG intentionally omitted. |
| Plotting/vignette/UI code | — | Intentionally omitted. |

## Dependency coverage

The root repository was checked before implementation. Compatible shared
packages are used through sibling FPM paths for `rfortran-core`,
`rfortran-linalg`, and `numDeriv`. Existing top-level translations also exist
for R dependencies such as `lme4`, `MASS`, and `Matrix`, but this package does
not duplicate them and does not add unused dependencies for R-only adapter code.
