# API map

`polyaAeppli` 2.0.2 exports four functions. All four are translated.

| R API | Fortran API |
|---|---|
| `dPolyaAeppli` | `d_polya_aeppli`, `d_polya_aeppli_vec` |
| `pPolyaAeppli` | `p_polya_aeppli`, `p_polya_aeppli_vec` |
| `qPolyaAeppli` | `q_polya_aeppli`, `q_polya_aeppli_vec` |
| `rPolyaAeppli` | `r_polya_aeppli`, `r_polya_aeppli_vec` |

The package's internal numerical helpers are also translated:

- `lPolyaAeppliArray` -> `log_pmf_array`
- `gArray` -> `log_cdf_array`
- `logTailPA` -> `log_tail_pa`
- `hArray` -> `log_sf_array`

The R-only helper `is.wholenumber` is folded into scalar density validation.
