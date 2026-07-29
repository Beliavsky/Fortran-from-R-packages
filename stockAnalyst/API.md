# API map

All functions are in module `stock_analyst`. `dp` is public.

| Original R export | Preferred Fortran procedure | Purpose |
|---|---|---|
| `shareValueUsingDDM1yr` | `share_value_using_ddm_1yr` | One-period DDM |
| `shareValueUsingDDMnYrs` | `share_value_using_ddm_n_years` | Multi-period DDM |
| `shareValueGGMconstantGrowth` | `share_value_ggm_constant_growth` | Constant-growth Gordon model |
| `shareValuePreferredStock` | `share_value_preferred_stock` | Perpetual preferred stock |
| `shareValueGGMNegativeGrowth` | `share_value_ggm_negative_growth` | Negative-growth Gordon model |
| `computingGusingGGM` | `computing_g_using_ggm` | Implied Gordon growth |
| `justifiedLeadingPE` | `justified_leading_pe` | Justified leading P/E |
| `justifiedTrailingPE` | `justified_trailing_pe` | Justified trailing P/E |
| `computingRwithGGM` | `computing_r_with_ggm` | Required return from Gordon model |
| `shareValUsingTwoStageDDM` | `share_val_using_two_stage_ddm` | Two-stage DDM terminal expression |
| `shareValUsingThreeStageDDM` | `share_val_using_three_stage_ddm` | Three-stage DDM terminal expression |
| `shareValUsingTwoStageHmodel` | `share_val_using_two_stage_hmodel` | Two-stage H-model |
| `shareValueNoCurrentDivdend` | `share_value_no_current_dividend` | Value before dividends begin |
| `annulizedHPR` | `annualized_hpr` | Annualized holding-period return |
| `firmValueUsingDiscFCFF` | `firm_value_using_disc_fcff` | Discounted FCFF firm value |
| `equityValueGivenDebtMV` | `equity_value_given_debt_mv` | FCFF equity value net of debt |
| `shareValueGivenDebtMV` | `share_value_given_debt_mv` | FCFF value per share |
| `shareValueUsingDiscFCFE` | `share_value_using_disc_fcfe` | Discounted FCFE per share |
| `firmValueConstantG` | `firm_value_constant_g` | Constant-growth FCFF firm value |
| `equityValueConstantG` | `equity_value_constant_g` | Constant-growth equity value |
| `shareValConstantG` | `share_val_constant_g` | Constant-growth FCFE per share |
| `shareValTwoStage` | `share_val_two_stage` | Two-stage FCFE vector valuation |
| `shareValThreeStg` | `share_val_three_stage` | Three-stage FCFE vector valuation |
| `shareValueRI` | `share_value_ri` | Residual-income value |
| `shareValueComputedRI` | `share_value_computed_ri` | Residual income from EPS/BVPS |
| `computingAbsRI` | `computing_abs_ri` | Absolute residual income from EBIT |
| `computingRI` | `computing_ri` | Residual-income vector |
| `shareValueROE` | `share_value_roe` | Residual-income value from ROE |
| `singleStageR` | `single_stage_r` | Single-stage residual-income value |
| `shareValueRImultiStageEPS` | `share_value_ri_multi_stage_eps` | Multistage RI from EPS |
| `shareValueRImultiStg` | `share_value_ri_multi_stage_roe` | Multistage RI from ROE |
| `shareValueRIplusPVTV` | `share_value_ri_plus_pvtv` | RI plus continuing value |
| `trailingPEbasicEPS` | `trailing_pe_basic_eps` | Trailing P/E using basic EPS |
| `trailingPEdilutedEPS` | `trailing_pe_diluted_eps` | Trailing P/E using diluted EPS |
| `earningYieldEP` | `earning_yield_ep` | Earnings yield |
| `leadingPEnext4Qs` | `leading_pe_next_4qs` | Leading P/E from four quarters |
| `leadingFY1PE` | `leading_fy1_pe` | FY1 leading P/E |
| `leadingFY2PE` | `leading_fy2_pe` | FY2 leading P/E |
| `predictedPEonCSR` | `predicted_pe_on_csr` | Cross-sectional regression P/E |
| `forwardPEG` | `forward_peg` | Forward PEG |
| `predictedPEbyFEDmodel` | `predicted_pe_by_fed_model` | Fed-model P/E |
| `impliedPEbyYardeniModel` | `implied_pe_by_yardeni_model` | Yardeni-model P/E |
| `sharePriceUsingPastPE` | `share_price_using_past_pe` | Share price from mean/median P/E |
| `PEforPassThroughInflation` | `pe_for_pass_through_inflation` | Inflation pass-through P/E |
| `terminalValueUsingPE` | `terminal_value_using_pe` | Comparable/GGM terminal value |
| `computingBVperShare` | `computing_bv_per_share` | Book value per common share |
| `computingPB` | `computing_pb` | Trailing/GGM price-to-book |
| `computingPS` | `computing_ps` | Trailing/GGM price-to-sales |
| `computingEVdollarVal` | `computing_ev_dollar_val` | Enterprise value |
| `computingEVmultiple` | `computing_ev_multiple` | EV/sales or EV/EBITDA |
| `computingSustainableG` | `computing_sustainable_g` | Sustainable growth |
| `computingRwithCAPM` | `computing_r_with_capm` | CAPM required return |
| `computingWACC` | `computing_wacc` | Weighted average cost of capital |
| `computingRwithFFM` | `computing_r_with_ffm` | Fama-French required return |
| `computingRwithHmodel` | `computing_r_with_hmodel` | H-model required return |

The original names are also exported as generic compatibility aliases. Fortran is case-insensitive, so the camel-case spelling shown above is optional.
