! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
program test_european
  use greeks
  implicit none
  type(greeks_result) :: call, put, cash, asset
  real(dp) :: h, up, down

  put=bs_european_greeks(120.0_dp,100.0_dp,0.02_dp,4.5_dp,0.22_dp,0.015_dp,payoff_put)
  call assert_close(put%fair_value,10.132463137873106_dp,1.0e-12_dp)
  call assert_close(put%delta,-0.23435480097375996_dp,1.0e-12_dp)
  call assert_close(put%gamma,0.005312008525061833_dp,1.0e-14_dp)
  call assert_close(put%vega,75.72799353328148_dp,1.0e-11_dp)

  call=bs_european_greeks(100.0_dp,105.0_dp,0.03_dp,1.5_dp,0.25_dp,0.01_dp,payoff_call)
  put=bs_european_greeks(100.0_dp,105.0_dp,0.03_dp,1.5_dp,0.25_dp,0.01_dp,payoff_put)
  call assert_close(call%fair_value-put%fair_value, &
    100.0_dp*exp(-0.01_dp*1.5_dp)-105.0_dp*exp(-0.03_dp*1.5_dp),1.0e-12_dp)

  cash=bs_european_greeks(100.0_dp,105.0_dp,0.03_dp,1.5_dp,0.25_dp,0.01_dp,payoff_cash_call)
  asset=bs_european_greeks(100.0_dp,105.0_dp,0.03_dp,1.5_dp,0.25_dp,0.01_dp,payoff_asset_call)
  call assert_close(asset%fair_value-105.0_dp*cash%fair_value,call%fair_value,1.0e-12_dp)

  h=1.0e-3_dp
  up=bs_european_price(100.0_dp+h,105.0_dp,0.03_dp,1.5_dp,0.25_dp,0.01_dp,payoff_cash_call)
  down=bs_european_price(100.0_dp-h,105.0_dp,0.03_dp,1.5_dp,0.25_dp,0.01_dp,payoff_cash_call)
  call assert_close(cash%delta,(up-down)/(2.0_dp*h),2.0e-10_dp)

  h=1.0e-4_dp
  up=bs_european_price(100.0_dp,105.0_dp,0.03_dp,1.5_dp,0.25_dp+h,0.01_dp,payoff_call)
  down=bs_european_price(100.0_dp,105.0_dp,0.03_dp,1.5_dp,0.25_dp-h,0.01_dp,payoff_call)
  call assert_close(call%vega,(up-down)/(2.0_dp*h),2.0e-6_dp)
  print '(a)', 'test_european: PASS'
contains
  subroutine assert_close(actual,expected,tol)
    real(dp),intent(in)::actual,expected,tol
    if (abs(actual-expected)>tol+tol*abs(expected)) then
      print '(a,3es25.16)', 'mismatch: ',actual,expected,abs(actual-expected)
      error stop 1
    end if
  end subroutine assert_close
end program test_european
