# API map

| R package routine | Modern Fortran equivalent |
|---|---|
| `CalcEAD` | `calc_ead` |
| `CalcPFE` | `calc_pfe` |
| `CalcRC` | `calc_rc` |
| `CalculateFactorMult` | `calculate_factor_multiplier` |
| `SingleTradeAddon` | `single_trade_addon` |
| `LoadSupervisoryData` | `load_supervisory_data`, `default_supervisory_data` |
| `HandleBasisVol` | `hedging_set_name` and the basis/volatility logic in `calc_addon` |
| `CreateTradeGraph` and `Group*Trades` | typed grouping performed by `calc_addon` |
| `CalcAddon` | `calc_addon` |
| `runExampleCalcs` | `calculate_portfolio` |
| `SACCRCalculator` | `saccr_calculator` |
| `DetermineCCRMethodology` | `determine_ccr_methodology` |
| `ExampleIRD` | `example_ird` |
| `ExampleFX` | `example_fx` |
| `ExampleCredit` | `example_credit` |
| `ExampleComm` | `example_commodity` |
| `ExampleIRDCredit` | `example_ird_credit` |
| `ExampleIRDCommMargined` | `example_ird_commodity_margined` |
| `ExampleBasisVol` | `example_basis_volatility` |
| `ExampleFXHedge` | `example_fx_hedge`, `apply_fx_hedge` |

## Result-tree mapping

The R `data.tree` hierarchy is represented as follows:

- portfolio root: `portfolio_result_t`
- counterparty/CSA netting set: `exposure_result_t`
- aggregate add-on node: `addon_result_t`
- asset-class nodes: `asset_class_result_t`
- hedging-set/reference-entity/currency nodes: `hedging_set_result_t`
- trade leaves and exposure details: `single_trade_addon_t`

This avoids runtime reflection and string-indexed tree mutation while retaining
the numerical hierarchy and allowing direct array processing.

## Intentionally omitted interfaces

- `data.tree` node construction and traversal
- JSON serialization through `jsonlite`
- R reference-class reflection
- package startup messages
- plotting or interactive display code

CSV input and exposure-summary CSV output are retained as native Fortran
routines.
