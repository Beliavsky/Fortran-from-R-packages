! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
program test_asian_implied
  use greeks
  implicit none
  type(greeks_result) :: g
  type(implied_vol_result) :: iv
  real(dp) :: target,h,up,down

  g=bs_geometric_asian_greeks(110.0_dp,100.0_dp,0.02_dp,4.5_dp,0.22_dp,0.015_dp,payoff_put)
  call assert_close(g%fair_value,6.566091953542262_dp,1.0e-12_dp)
  call assert_close(g%delta,-0.2920249031792912_dp,1.0e-12_dp)
  call assert_close(g%gamma,0.01097609568886642_dp,1.0e-13_dp)

  h=1.0e-4_dp
  up=bs_geometric_asian_price(110.0_dp,100.0_dp,0.02_dp,4.5_dp,0.22_dp+h,0.015_dp,payoff_put)
  down=bs_geometric_asian_price(110.0_dp,100.0_dp,0.02_dp,4.5_dp,0.22_dp-h,0.015_dp,payoff_put)
  call assert_close(g%vega,(up-down)/(2.0_dp*h),2.0e-5_dp)

  target=bs_european_price(100.0_dp,100.0_dp,0.03_dp,5.0_dp,0.27_dp,0.015_dp,payoff_call)
  iv=bs_implied_volatility(target,100.0_dp,100.0_dp,0.03_dp,5.0_dp,0.015_dp,payoff_call)
  if (.not.iv%converged) error stop 1
  call assert_close(iv%volatility,0.27_dp,1.0e-12_dp)

  target=bs_geometric_asian_price(100.0_dp,95.0_dp,0.02_dp,2.0_dp,0.31_dp,0.01_dp,payoff_call)
  iv=geometric_asian_implied_volatility(target,100.0_dp,95.0_dp,0.02_dp,2.0_dp,0.01_dp,payoff_call)
  if (.not.iv%converged) error stop 1
  call assert_close(iv%volatility,0.31_dp,2.0e-7_dp)
  print '(a)', 'test_asian_implied: PASS'
contains
  subroutine assert_close(actual,expected,tol)
    real(dp),intent(in)::actual,expected,tol
    if (abs(actual-expected)>tol+tol*abs(expected)) then
      print '(a,3es25.16)', 'mismatch: ',actual,expected,abs(actual-expected)
      error stop 1
    end if
  end subroutine assert_close
end program test_asian_implied
