# R-to-Fortran API map

| R export | Fortran routine | Translation status |
|---|---|---|
| `bond_future_price` | `yrnd_api::bond_future_price` | Computational output translated; plots and descriptive labels omitted. |
| `stir_future_price` | `yrnd_api::stir_future_price` | Computational output translated; plots and labels omitted. |
| `stir_rate` | `yrnd_api::stir_rate` | Price-to-rate density transformation translated. |
| `ctd_bond_yield` | `yrnd_api::ctd_bond_yield` | Bond cash flows, carry, `xirr`, Jacobian density transformation translated. |
| `proba_ctd` | `yrnd_api::proba_ctd` | CTD self-consistency and probability integration translated. |
| `proba_ctd_opt` | `yrnd_api::proba_ctd_opt` | Option-to-futures carry and CTD probabilities translated. |
| `bond_yield_spread` | `yrnd_api::bond_yield_spread` | Gaussian-copula simulation, CTD yields, spread samples, and KDE translated. |
| `bond_future_charac_bbg` | none | Bloomberg data retrieval omitted. |
| `stir_future_charac_bbg` | none | Bloomberg data retrieval omitted. |
| `deliv_bonds_charac_bbg` | none | Bloomberg data retrieval omitted. |
| `option_prices_bbg` | none | Bloomberg data retrieval omitted. |

## Result types

- `density_result_t` replaces the R list returned for future-price densities.
- `transformed_density_t` replaces the R lists returned for rate/yield densities.
- `ctd_probability_result_t` replaces the CTD probability data frame.
- `spread_result_t` contains simulated spread observations and a Gaussian KDE.
- `bond_t`, `bond_context_t`, and `date_t` replace heterogeneous data frames.
