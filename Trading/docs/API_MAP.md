# R to Fortran API map

## Dependence and entropy

| R API | Fortran API |
|---|---|
| `AngularDistance` | `angular_distance` |
| `Chebyshev_distance` | `chebyshev_distance` |
| `SampleEntropy` | `sample_entropy` |
| `CrossSampleEntropy` | `cross_sample_entropy` |
| `NormXASampEn` | `normalized_cross_sample_entropy` |
| `VariationOfInformation` | `variation_of_information` |
| `InformationAdjustedCorr` | `information_adjusted_correlation` |
| `InformationAdjustedBeta` | `information_adjusted_beta` |
| `DynamicBeta` | `dynamic_beta` and `dynamic_beta_result` |

## Trade and exposure objects

The R reference-class hierarchy is represented by one tagged derived type,
`trade_t`. Call `trade%configure_class("IRDSwap")`, for example, then populate
its public fields.

| R class or method | Fortran equivalent |
|---|---|
| `Trade`, `IRD`, `FX`, `Credit`, `Commodity`, `Equity`, subclasses | `trade_t` plus `configure_class` |
| `CalcAdjNotional` | `trade%calc_adjusted_notional()` |
| `CalcSupervDuration` | `trade%calc_supervisory_duration()` |
| `CalcMaturityFactor` | `trade%calc_maturity_factor()` |
| `CalcSupervDelta` | `trade%calc_supervisory_delta()` |
| `SetTimeBucket` | `trade%set_time_bucket()` |
| `isBasisSwap` | `trade%is_basis_swap()` |
| `ComputeVarianceUnits` | `trade%compute_variance_units()` |
| `CalcMaturity` | `trade%calc_option_maturity()` |
| `setFXDynamic` | `trade%set_fx_dynamic()` |
| `SplitBondExposure` | `split_bond_exposure` |
| `SelectDerivatives` | `select_derivatives` |
| `ParseTrades` | `parse_trades_csv` |
| `GetTradeDetails` | `write_trade_details` or direct field access |
| `CSA` | `csa_t` |
| `CSA$ApplyThres` | `csa%apply_threshold()` |
| `CSA$CalcMF` | `csa%maturity_factor()` |
| `Collateral` | `collateral_t` |
| `Curve` | `curve_t` |
| `Curve$CalcInterpPoints` | `curve%interpolate()` |
| `HashTable` | `hash_table_t` |

## Climate metrics

| R API | Fortran API |
|---|---|
| `Carbon_Footprint` | `carbon_footprint` |
| `Carbon_Intensity` | `carbon_intensity` |
| `Total_Carbon_Emissions` | `total_carbon_emissions` |
| `Weighted_Average_Carbon_Intensity` | `weighted_average_carbon_intensity` |

The Fortran routines take issuer-aligned numeric arrays. Any issuer-name join is
performed by the caller before invocation.

## Lotteries

| R API | Fortran API |
|---|---|
| `EuroMillionsResults`, `EuroJackpotResults` | `load_euro_lottery_results` |
| `SetForLifeResults` | `load_set_for_life_results` |
| `UKLotteryResults` | `load_uk_lottery_results` |
| `UKThunderBallResults` | `load_uk_thunderball_results` |
| `EuroLotteryBacktesting` | `euro_lottery_backtest` |
| `SetForLifeBacktesting` | `set_for_life_backtest` |
| `UKLotteryBacktesting` | `uk_lottery_backtest` |
| `UKThunderballBacktesting` | `uk_thunderball_backtest` |
| `CalcEuroLotteryPnL` | `calculate_euro_lottery_pnl` |
| `CalcSetForLifePnL` | `calculate_set_for_life_pnl` |
| `CalcUKLotteryPnL` | `calculate_uk_lottery_pnl` |
| `CalcUKThunderBallPnL` | `calculate_uk_thunderball_pnl` |
| `top5` | `top_five_numbers` |
| `EuroLotteryAllCombinations` | `euro_combination_iterator_t` |
| `OuterJoinMerge` | `outer_join_merge_integer` |

## Betting simulations

| R API | Fortran API |
|---|---|
| `capped_fibonacci_seq` | `capped_fibonacci_sequence` |
| `martingale_strategy_repetitions` | `martingale_strategy_repetitions` |
| `roulette_pl_calculator_dalembert` | `roulette_dalembert` |
| `roulette_pl_calculator_fibonacci` | `roulette_fibonacci` |
| `roulette_pl_calculator_labouchere` | `roulette_labouchere` |
| `roulette_pl_calculator_martingale` | `roulette_martingale` |
| `roulette_pl_calculator_specific_number` | `roulette_specific_number` |

All simulation routines accept an optional integer seed and return arrays in a
`betting_result_t` or `repetitions_result_t` object. Plotting is omitted.
