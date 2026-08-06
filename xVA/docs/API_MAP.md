# R to Fortran API map

| R function | Modern Fortran routine/type | Module |
|---|---|---|
| `CalcNGR` | `calc_ngr` | `xva_core` |
| `CalcPD` | `calc_pd` | `xva_core` |
| `CalcVA` | `calc_va` | `xva_core` |
| `GenerateTimeGrid` | `generate_time_grid` | `xva_core` |
| `IS_ELIGIBLE_CCY` | `is_eligible_currency` | `xva_core` |
| `IS_IG` | `is_investment_grade` | `xva_core` |
| `CalcSimulatedExposure` | `calc_simulated_exposure` | `xva_exposure` |
| `calcEADRegulatory` | `calc_ead_regulatory` | `xva_regulatory` |
| `calcEffectiveMaturity` | `calc_effective_maturity` | `xva_core` |
| `calcDefCapital` | `calc_default_capital` | `xva_core` |
| `calcCVACapital` | `calc_cva_capital` | `xva_regulatory` |
| `calcKVA` | `calc_kva` | `xva_core` |
| `LoadSupervisoryCVAData` | `load_supervisory_cva_data` | `xva_supervisory` |
| `xVACalculator` | `xva_calculator` | `xva_calculator_mod` |
| `xVACalculatorExample` | `app/xva_demo.f90` | example program |
| R simulation list | `simulation_data_t` | `xva_types` |
| R regulatory list | `regulatory_data_t` | `xva_types` |
| R exposure-profile list | `exposure_profile_t` | `xva_types` |
| R CVA-sensitivity list | `cva_sensitivity_t` | `xva_types` |
| R xVA result list | `xva_result_t` | `xva_types` |

The convenience module `xva` re-exports the public API and the required
`Trading` types `trade_t`, `csa_t`, `collateral_t`, and `curve_t`.
