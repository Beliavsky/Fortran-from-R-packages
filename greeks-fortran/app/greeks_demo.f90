! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
program greeks_demo
  use greeks
  implicit none
  type(greeks_result) :: euro, american, asian
  type(implied_vol_result) :: iv
  real(dp) :: target

  euro=bs_european_greeks(100.0_dp,105.0_dp,0.03_dp,1.5_dp,0.25_dp,0.01_dp,payoff_call)
  american=binomial_american_greeks(100.0_dp,105.0_dp,0.03_dp,1.5_dp,0.25_dp,0.01_dp,payoff_put,500)
  asian=bs_geometric_asian_greeks(100.0_dp,105.0_dp,0.03_dp,1.5_dp,0.25_dp,0.01_dp,payoff_call)
  target=bs_european_price(100.0_dp,105.0_dp,0.03_dp,1.5_dp,0.25_dp,0.01_dp,payoff_call)
  iv=bs_implied_volatility(target,100.0_dp,105.0_dp,0.03_dp,1.5_dp,0.01_dp,payoff_call)

  print '(a,f12.6)', 'European call value:      ',euro%fair_value
  print '(a,f12.6)', 'European call delta:      ',euro%delta
  print '(a,f12.6)', 'European call gamma:      ',euro%gamma
  print '(a,f12.6)', 'American put value:       ',american%fair_value
  print '(a,f12.6)', 'Geometric Asian call:     ',asian%fair_value
  print '(a,f12.6)', 'Recovered volatility:     ',iv%volatility
end program greeks_demo
