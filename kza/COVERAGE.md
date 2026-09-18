# API coverage

The translation status is **substantial**, with **9 of 9 (100%)** exported computational R functions mapped to public Fortran procedures. The percentage is function-mapping coverage, not a claim of complete R-interface compatibility.

| R function | Status | Main compatibility notes |
|---|---|---|
| `kz` | substantial | Numerical 1D/2D/3D iterated moving averages translated; R `ts` attributes are omitted. |
| `kza` | substantial | Adaptive 1D/2D/3D kernels, normalization modes, matrix symmetrization, and 1D tail marking translated; S3 list/call construction is omitted. |
| `kzsv` | substantial | Adaptive local sample variance translated; numeric inputs replace the R `kza` object. |
| `kzs` | substantial | Zero-frequency KZFT and automatic window rule translated; R `ts` metadata omitted. |
| `transfer_function` | complete | Direct numerical translation including the default frequency grid. |
| `kzft` | substantial | Missing-aware iterative local transform translated, including even/non-integer window alignment; R attributes omitted. |
| `kztp` | substantial | Frequency bank and third-order averaging translated; progress printing omitted. |
| `periodogram` | substantial | Same transform convention and returned magnitudes; direct DFT replaces FFT for dependency-free portability. |
| `rlv` | complete | 1D/2D/3D clamp and zero padding plus the `krnl=1` legacy corner-window definition translated. |

## Excluded from the denominator

Plotting functions and presentation-only registered methods are excluded. `summary.kzp` and `summary.kzsv` primarily select and print periods/dates of interest and are treated as reporting behavior rather than exported numerical kernels. Internal unexported helpers are likewise outside the exported-function denominator.

## Deliberate interface differences

Fortran returns numerical arrays instead of R S3 lists and does not reproduce R calls, `ts` attributes, plotting hooks, or console messages. `kz`, `kza`, and `rlv` use generic rank-specific procedures. Anisotropic KZ/KZA windows use optional `m2` and `m3` scalar arguments rather than an R integer vector.
