# Translation coverage

| R export | Fortran interface | Status |
|---|---|---|
| `g_ret` | generic `g_ret` | Implemented |
| `tail_mixture` | `tail_mixture` | Implemented |
| `infoset` | `infoset_estimate` | Implemented |
| `create_overlapping_windows` | `create_overlapping_windows` | Implemented |
| `LR_cp` | `lr_cp` | Implemented |
| `ptf_construction` | `ptf_construction` | Implemented |
| `summary_ptf` | `summary_ptf` | Implemented |
| `plot_ptf` | none | Plotting omitted |

The non-exported `plot_LR_cp` helper and all graphics calls inside
`tail_mixture` are omitted. Bundled R `.rda` datasets and the rendered vignette
are not included in the compiled project.

Dependencies replaced internally:

- `mixtools::normalmixEM`: deterministic two-normal EM.
- `quadprog::solve.QP`: modernized Goldfarb-Idnani solver.
- `Matrix::nearPD`: eigenvalue-clipped nearest positive-definite matrix.
- `stats::density`, `hist`, and `quantile`: explicit numerical routines.
