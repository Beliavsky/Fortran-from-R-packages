# API coverage

This file complements the machine-readable `[extra.translation]` tables in `fpm.toml`.

## Exported computational API

The coverage denominator is the 11 exported computational functions/methods identified from the supplied optR 1.2.5 `NAMESPACE`. `print.optR` and `summary.optR` are excluded as presentation-only methods.

All 11 have meaningful Fortran mappings. `optR.formula` is partial because Fortran receives an already constructed numeric design matrix rather than parsing R formulas and model-frame metadata. Other mapped S3 entry points are substantial rather than complete when their numerical operation is present but R class/list construction is not.

## Translated internal numerical kernels

The Fortran implementation also covers numerical routines that are not part of the exported-function denominator:

| Upstream R routine | Fortran procedure or implementation |
|---|---|
| `LU.optR` | `lu_solve`, `lu_decompose` |
| `LU.decompose` | `lu_decompose` |
| `cgm` | `conjugate_gradient` |
| `choleskilm` | `cholesky_solve` |
| `forwardsubsitution.Cholesky` | internal forward solve in `cholesky_solve` |
| `optR.gauss` | `gaussian_eliminate`, `gaussian_solve` |
| `opt.matrix.reorder` | `matrix_reorder` |
| `optR.multiplyfactor` | `multiply_factor` |
| `gaussSeidel` | `gauss_seidel` |
| `nonDiagMultipication` | `non_diag_multiplication` |
| `machinePrecision` | `zero_small_values` |
| `newtonRapson` | `newton_raphson` |
| `optR.fit` | `optr_solve`, `optr_fit_design` |
| `updateWeightsobs` | weight transformation inside `optr_fit_design` |
| `optR.backsubsitution` | `back_substitution` |
| `forwardsubsitution.optR` | `forward_substitution` |

The small `assign.values` helper only rearranges R list members and has no independent numerical operation.

## Material differences from R

- Formula parsing, `model.frame`, `model.matrix`, factor levels, contrasts, NA actions, R recycling/coercion, and S3 class metadata are not reproduced.
- Plotting and printing are skipped.
- Errors and warnings are represented by integer status codes.
- LU/Gaussian pivot handling is made deterministic and numerically safer; this can differ from upstream behavior for singular or nearly singular matrices.
- The Jacobian keeps the upstream forward-difference default step of `1e-4` and assumes a square nonlinear system, as the R implementation does.
- `optr_fit_design` follows upstream rectangular fitting by forming normal equations. It does not silently switch to QR/SVD, because that would change the algorithm being translated.
