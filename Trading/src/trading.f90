module trading
  use trading_kinds, only : dp, str_len
  use trading_stats, only : mean_dp, variance_dp, sd_dp, covariance_dp, &
    correlation_dp, correlation_matrix, quantile_type7, normal_cdf, &
    set_random_seed, random_normal
  use trading_dependence, only : angular_distance, chebyshev_distance, &
    sample_entropy, cross_sample_entropy, normalized_cross_sample_entropy, &
    variation_of_information, information_adjusted_correlation, &
    information_adjusted_beta
  use trading_dynamic_beta, only : dynamic_beta_result, dynamic_beta
  use trading_trades, only : trade_t, split_bond_exposure, select_derivatives
  use trading_csa, only : csa_t, collateral_t
  use trading_curve, only : curve_t, natural_cubic_spline, linear_interpolation
  use trading_hash_table, only : hash_table_t
  use trading_climate, only : carbon_footprint, carbon_intensity, &
    total_carbon_emissions, weighted_average_carbon_intensity
  use trading_lottery, only : lottery_draw_t, lottery_pnl_result_t, &
    euro_combination_iterator_t, load_euro_lottery_results, &
    load_set_for_life_results, load_uk_lottery_results, &
    load_uk_thunderball_results, euro_lottery_backtest, &
    set_for_life_backtest, uk_lottery_backtest, uk_thunderball_backtest, &
    calculate_euro_lottery_pnl, calculate_set_for_life_pnl, &
    calculate_uk_lottery_pnl, calculate_uk_thunderball_pnl, &
    top_five_numbers, euro_lottery_combination_count, &
    outer_join_merge_integer, parse_date
  use trading_betting, only : betting_result_t, repetitions_result_t, &
    capped_fibonacci_sequence, martingale_strategy_repetitions, &
    roulette_dalembert, roulette_fibonacci, roulette_labouchere, &
    roulette_martingale, roulette_specific_number
  use trading_io, only : parse_trades_csv, load_curve_csv, load_csa_csv, &
    load_collateral_csv, load_track_record_csv, write_trade_details
  implicit none
  public

end module trading
