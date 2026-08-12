module matchingr
   use matchingr_kinds, only : dp
   use matchingr_types
   use matchingr_utils, only : sort_index, sort_index_one_sided, rank_index, &
      check_preferences, check_roommate_preferences, repeat_rows_real, repeat_cols_real
   use matchingr_galeshapley, only : marriage_market, marriage_market_preferences, &
      college_admissions, college_admissions_preferences, gale_shapley_stable
   use matchingr_roommate, only : stable_roommates, stable_roommates_preferences, roommate_stable
   use matchingr_ttc, only : top_trading_cycles, top_trading_cycles_preferences, top_trading_stable
   implicit none
   public
end module matchingr
