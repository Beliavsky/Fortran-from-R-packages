# API map

| R export | Fortran API | Status |
|---|---|---|
| `TP` | `tp`, `touchard_polynomial` | translated |
| `bell_mle` | `bell_mle` | translated |
| `mle.bell` | `mle_bell` | translated |
| `mle_borel` | `mle_borel` | translated |
| `mle_poisson` | `mle_poisson` | translated |
| `dbellt/pbellt/qbellt/rbellt` | same names | translated |
| `mle_bt` | `mle_bt` | translated |
| `dzibellt/pzibellt/qzibellt/rzibellt` | same names | translated |
| `mle_zibellt` | `mle_zibellt` | translated |
| `mle_zibell` | `mle_zibell` | translated |
| `mle_zip` | `mle_zip` | translated |
| `mle_zoibell` | `mle_zoibell` | translated |
| `mle_zoip` | `mle_zoip` | translated |
| `data_criminal` | `data_criminal()` | translated |
| `data_sbirth` | `data_sbirth()` | translated |

The R package contains no plotting routines or compiled C/C++/Fortran code.
Its R list/column-name presentation layer is replaced by `mle_result_t`.
