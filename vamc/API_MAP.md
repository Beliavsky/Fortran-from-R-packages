# API map

| Upstream R routine | Fortran routine/type | Status |
|---|---|---|
| `buildCurve` | `build_curve`, `build_curve_named`, `yield_curve_type` | Implemented |
| `genIndexScen` | `gen_index_scen` | Implemented |
| `genFundScen` | generic `gen_fund_scen` for rank-2/rank-3 arrays | Implemented |
| `genPortInception` | `gen_port_inception`, `policy_type`, `portfolio_type` | Implemented |
| `calcMortFactors` | `calc_mort_factors`, `mortality_factors_type` | Implemented |
| `valuateOnePolicy` | generic `valuate_one_policy` | Implemented |
| `valuatePortfolio` | `valuate_portfolio` | Implemented |
| `ageOnePolicy` | `age_one_policy` | Implemented |
| `agePortfolio` | `age_portfolio` | Implemented |
| private `projectDBRP` ... `projectDBWB` | `project_policy` dispatch on `policy_type%product_type` | All 19 implemented |
| private date/calendar helpers | `vamc_dates` module | Implemented |
| private swap helpers | `log_linear_discount`, `present_value_swap` | Implemented |
| private `rFundMap` | `random_fund_map` | Implemented |
| `readMortTable` | `make_mortality_table`; file parsing left to caller | Numerical content implemented |
| package data objects | Original `.rda` files under `upstream/vamc-master/data` | Preserved, not converted |
