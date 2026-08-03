# API reference

All public procedures are available through:

```fortran
use fincal
```

## R-to-Fortran name mapping

| R name | Fortran name |
|---|---|
| `EIR` | `eir` |
| `EPS` | `eps` |
| `SFRatio` | `sfratio` or `sf_ratio` |
| `Sharpe.ratio` | `sharpe_ratio` |
| `bdy` | `bdy` |
| `bdy2mmy` | `bdy2mmy` or `bdy_to_mmy` |
| `cash.ratio` | `cash_ratio` |
| `coefficient.variation` | `coefficient_variation` |
| `cogs` | `cogs` |
| `current.ratio` | `current_ratio` |
| `ddb` | `double_declining_balance` |
| `debt.ratio` | `debt_ratio` |
| `diluted.EPS` | `diluted_eps` |
| `discount.rate` | `discount_rate` / `solve_discount_rate` |
| `ear` | `ear` |
| `ear.continuous` | `ear_continuous` |
| `ear2bey` | `ear2bey` or `ear_to_bey` |
| `ear2hpr` | `ear2hpr` or `ear_to_hpr` |
| `financial.leverage` | `financial_leverage` |
| `fv` | `fv` |
| `fv.annuity` | `fv_annuity` |
| `fv.simple` | `fv_simple` |
| `fv.uneven` | `fv_uneven` |
| `geometric.mean` | `geometric_mean` |
| `gpm` | `gpm` or `gross_profit_margin` |
| `harmonic.mean` | `harmonic_mean` |
| `hpr` | `hpr` |
| `hpr2bey` | `hpr2bey` or `hpr_to_bey` |
| `hpr2ear` | `hpr2ear` or `hpr_to_ear` |
| `hpr2mmy` | `hpr2mmy` or `hpr_to_mmy` |
| `irr` | `irr` / `solve_irr` |
| `irr2` | `irr2` |
| `iss` | `iss` or `issuable_shares` |
| `lt.d2e` | `lt_d2e` or `long_term_debt_to_equity` |
| `mmy2hpr` | `mmy2hpr` or `mmy_to_hpr` |
| `n.period` | `n_period` |
| `npm` | `npm` or `net_profit_margin` |
| `npv` | `npv` |
| `pmt` | `pmt` |
| `pv` | `pv` |
| `pv.annuity` | `pv_annuity` |
| `pv.perpetuity` | `pv_perpetuity` |
| `pv.simple` | `pv_simple` |
| `pv.uneven` | `pv_uneven` |
| `quick.ratio` | `quick_ratio` |
| `r.continuous` | `continuous_rate` or `r_continuous` |
| `r.norminal` | `nominal_rate` or compatibility alias `r_norminal` |
| `r.perpetuity` | `r_perpetuity` or `perpetuity_rate` |
| `sampling.error` | `sampling_error` |
| `slde` | `slde` or `straight_line_depreciation` |
| `total.d2e` | `total_d2e` or `total_debt_to_equity` |
| `twrr` | `twrr` |
| `was` | `was` or `weighted_average_shares` |
| `wpr` | `wpr` or `weighted_portfolio_return` |

## Result types and status codes

`cogs` returns `type(inventory_result)` with:

- `cost_of_goods`
- `ending_inventory`
- `status`

`solve_discount_rate` and `solve_irr` return `type(root_result)` with:

- `root`
- `function_value`
- `iterations`
- `status`

Status constants are `fincal_ok`, `fincal_invalid_input`, `fincal_size_mismatch`, `fincal_no_root`, `fincal_insufficient_inventory`, `fincal_nonfinite_result`, and `fincal_weights_not_unit`.

## Sign convention

The translation preserves FinCal's TVM sign convention: cash paid and cash received have opposite signs. Thus a positive future value generally corresponds to a negative present value or payment.
