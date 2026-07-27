! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
program test_malliavin
  use greeks
  implicit none
  type(mc_greeks_result) :: mc, geom, asian, cv
  type(greeks_result) :: exact
  real(dp) :: allowance

  mc=malliavin_european_greeks(100.0_dp,100.0_dp,0.03_dp,1.0_dp,0.2_dp, &
    payoff_call,40000,123,.true.)
  exact=bs_european_greeks(100.0_dp,100.0_dp,0.03_dp,1.0_dp,0.2_dp,0.0_dp,payoff_call)
  allowance=7.0_dp*mc%standard_error%fair_value+0.02_dp
  if (abs(mc%estimate%fair_value-exact%fair_value)>allowance) error stop 1
  allowance=7.0_dp*mc%standard_error%delta+0.005_dp
  if (abs(mc%estimate%delta-exact%delta)>allowance) error stop 1

  geom=malliavin_geometric_asian_greeks(100.0_dp,100.0_dp,0.02_dp,1.0_dp,0.25_dp, &
    0.0_dp,payoff_call,32,12000,77,.true.)
  exact=bs_geometric_asian_greeks(100.0_dp,100.0_dp,0.02_dp,1.0_dp,0.25_dp,0.0_dp,payoff_call)
  allowance=8.0_dp*geom%standard_error%fair_value+0.05_dp
  if (abs(geom%estimate%fair_value-exact%fair_value)>allowance) error stop 1

  asian=malliavin_asian_greeks(100.0_dp,100.0_dp,0.02_dp,1.0_dp,0.25_dp, &
    0.0_dp,payoff_call,32,10000,91,.false.)
  cv=bs_malliavin_asian_greeks(100.0_dp,100.0_dp,0.02_dp,1.0_dp,0.25_dp, &
    0.0_dp,payoff_call,32,10000,91)
  if (.not.asian%ok .or. .not.cv%ok) error stop 1
  if (asian%estimate%fair_value<=0.0_dp .or. cv%estimate%fair_value<=0.0_dp) error stop 1
  if (cv%standard_error%fair_value>=asian%standard_error%fair_value) error stop 1
  if (abs(asian%estimate%delta_d-asian%estimate%delta)>0.25_dp) error stop 1
  print '(a)', 'test_malliavin: PASS'
end program test_malliavin
