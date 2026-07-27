! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
program analytic_greeks
  use greeks
  implicit none
  type(greeks_result) :: result

  result=bs_european_greeks(120.0_dp,100.0_dp,0.02_dp,4.5_dp,0.22_dp,0.015_dp,payoff_put)
  if (.not.result%ok) error stop trim(result%message)
  print '(a,f12.6)', 'fair value = ',result%fair_value
  print '(a,f12.6)', 'delta      = ',result%delta
  print '(a,f12.6)', 'vega       = ',result%vega
  print '(a,f12.6)', 'gamma      = ',result%gamma
  print '(a,f12.6)', 'vomma      = ',result%vomma
end program analytic_greeks
