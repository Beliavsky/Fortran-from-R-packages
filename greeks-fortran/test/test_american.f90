! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
program test_american
  use greeks
  implicit none
  type(greeks_result) :: g
  type(implied_vol_result) :: iv
  real(dp) :: european,target

  g=binomial_american_greeks(100.0_dp,100.0_dp,0.05_dp,1.0_dp,0.2_dp,0.0_dp,payoff_put,500)
  call assert_close(g%fair_value,6.092808546444233_dp,3.0e-10_dp)
  if (g%delta>=0.0_dp .or. g%gamma<=0.0_dp .or. g%vega<=0.0_dp) error stop 1
  european=bs_european_price(100.0_dp,100.0_dp,0.05_dp,1.0_dp,0.2_dp,0.0_dp,payoff_put)
  if (g%fair_value<european) error stop 1

  g=binomial_american_greeks(100.0_dp,100.0_dp,0.05_dp,1.0_dp,0.2_dp,0.0_dp,payoff_call,300)
  european=bs_european_price(100.0_dp,100.0_dp,0.05_dp,1.0_dp,0.2_dp,0.0_dp,payoff_call)
  call assert_close(g%fair_value,european,2.0e-10_dp)

  target=binomial_american_price(100.0_dp,105.0_dp,0.03_dp,1.25_dp,0.28_dp,0.01_dp,payoff_put,350)
  iv=american_implied_volatility(target,100.0_dp,105.0_dp,0.03_dp,1.25_dp,0.01_dp,payoff_put,350)
  if (.not.iv%converged) error stop 1
  call assert_close(iv%volatility,0.28_dp,2.0e-6_dp)
  print '(a)', 'test_american: PASS'
contains
  subroutine assert_close(actual,expected,tol)
    real(dp),intent(in)::actual,expected,tol
    if (abs(actual-expected)>tol+tol*abs(expected)) then
      print '(a,3es25.16)', 'mismatch: ',actual,expected,abs(actual-expected)
      error stop 1
    end if
  end subroutine assert_close
end program test_american
