! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
program malliavin_greeks
  use greeks
  implicit none
  type(mc_greeks_result) :: result

  result=bs_malliavin_asian_greeks(100.0_dp,100.0_dp,0.02_dp,1.0_dp,0.25_dp, &
    0.0_dp,payoff_call,32,20000,42)
  if (.not.result%ok) error stop trim(result%message)
  print '(a,f12.6,a,f10.6)', 'Asian call value = ',result%estimate%fair_value, &
    '  SE = ',result%standard_error%fair_value
  print '(a,f12.6,a,f10.6)', 'Asian call delta = ',result%estimate%delta, &
    '  SE = ',result%standard_error%delta
  print '(a,f12.6,a,f10.6)', 'Asian call vega  = ',result%estimate%vega, &
    '  SE = ',result%standard_error%vega
end program malliavin_greeks
