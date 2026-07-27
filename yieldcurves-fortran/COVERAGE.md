# Computational coverage

Upstream package: `yieldcurves` 0.1.0.

| Upstream routine | Fortran routine | Status |
|---|---|---|
| `yc_curve` | `yc_curve` | Complete |
| `yc_nelson_siegel` | `yc_nelson_siegel` | Complete |
| `yc_svensson` | `yc_svensson` | Complete |
| `yc_cubic_spline` | `yc_cubic_spline` | Natural and FMM |
| `yc_fit` | `yc_fit` | Complete dispatch |
| `yc_predict` | `yc_predict` | Complete |
| `yc_interpolate` | `yc_interpolate` | Linear, log-linear, cubic |
| `yc_discount` | `yc_discount` | Continuous, annual, semi-annual |
| `yc_forward` | `yc_forward` | Instantaneous and forward-forward |
| `yc_par_to_zero` | `yc_par_to_zero` | Annual and semi-annual |
| `yc_zero_to_par` | `yc_zero_to_par` | Annual and semi-annual |
| `yc_duration` | `yc_duration` | Complete |
| `yc_bond_duration` | `yc_bond_duration` | Complete |
| `yc_zspread` | `yc_zspread` | Complete |
| `yc_key_rate_duration` | `yc_key_rate_duration` | Complete |
| `yc_carry` | `yc_carry` | Complete |
| `yc_slope` | `yc_slope` | Complete; unavailable tenors return NaN |
| `yc_level_slope_curvature` | `yc_level_slope_curvature` | Complete |
| `yc_pca` | `yc_pca` | Complete numerical output |

The analytical helper formulas for Nelson-Siegel and Svensson rates, loadings,
and instantaneous forward rates are also public.

## Excluded presentation infrastructure

- `print.yc_curve`, `summary.yc_curve`, and `plot.yc_curve`
- `print.yc_pca`, `summary.yc_pca`, and `plot.yc_pca`
- `cli`, `graphics`, R dates, data frames, and S3 class machinery

The arrays needed for reporting and plotting are returned through typed
Fortran result structures.
