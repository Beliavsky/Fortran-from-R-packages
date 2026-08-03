! SPDX-License-Identifier: GPL-2.0-or-later
program test_rates_ratios
   use fincal
   implicit none
   real(dp), parameter :: tol = 2.0e-12_dp

   call assert_close(eir(0.05_dp, 1, 12), 0.00407412378365_dp, tol, 'eir')
   call assert_close(eir(0.05_dp, p = 12, rate_type = 'p'), 0.00416666666667_dp, tol, 'eir proportional')
   call assert_close(ear(0.12_dp, 12), 0.126825030132_dp, tol, 'ear')
   call assert_close(ear_continuous(0.1_dp), 0.105170918075648_dp, tol, 'ear continuous')
   call assert_close(ear_to_bey(0.08_dp), 0.078460969083_dp, tol, 'ear_to_bey')
   call assert_close(ear_to_hpr(0.05039_dp, 150.0_dp), 0.020408835325514_dp, 1.0e-9_dp, 'ear_to_hpr')
   call assert_close(bdy(1500.0_dp, 100000.0_dp, 120.0_dp), 0.045_dp, tol, 'bdy')
   call assert_close(bdy_to_mmy(0.045_dp, 120.0_dp), 0.045685279188_dp, tol, 'bdy_to_mmy')
   call assert_close(hpr(33.0_dp, 30.0_dp, 0.5_dp), 0.116666666667_dp, tol, 'hpr')
   call assert_close(hpr_to_bey(0.02_dp, 3.0_dp), 0.0808_dp, tol, 'hpr_to_bey')
   call assert_close(hpr_to_ear(0.015228_dp, 120.0_dp), 0.047042340436985_dp, tol, 'hpr_to_ear')
   call assert_close(hpr_to_mmy(0.01523_dp, 120.0_dp), 0.04569_dp, tol, 'hpr_to_mmy')
   call assert_close(mmy_to_hpr(0.04898_dp, 150.0_dp), 0.020408333333_dp, tol, 'mmy_to_hpr')
   call assert_close(nominal_rate(continuous_rate(0.03_dp, 4), 4), 0.03_dp, tol, 'continuous round trip')

   call assert_close(sf_ratio(0.09_dp, 0.03_dp, 0.12_dp), 0.5_dp, tol, 'sf_ratio')
   call assert_close(sharpe_ratio(0.038_dp, 0.015_dp, 0.07_dp), 0.328571428571_dp, tol, 'sharpe')
   call assert_close(coefficient_variation(0.15_dp, 0.39_dp), 0.384615384615385_dp, tol, &
      'coefficient variation')
   call assert_close(cash_ratio(3000.0_dp, 2000.0_dp, 2000.0_dp), 2.5_dp, tol, 'cash_ratio')
   call assert_close(current_ratio(8000.0_dp, 2000.0_dp), 4.0_dp, tol, 'current_ratio')
   call assert_close(quick_ratio(3000.0_dp, 2000.0_dp, 1000.0_dp, 2000.0_dp), 3.0_dp, tol, 'quick_ratio')
   call assert_close(debt_ratio(6000.0_dp, 20000.0_dp), 0.3_dp, tol, 'debt_ratio')
   call assert_close(long_term_debt_to_equity(8000.0_dp, 20000.0_dp), 0.4_dp, tol, 'lt debt/equity')
   call assert_close(total_debt_to_equity(6000.0_dp, 20000.0_dp), 0.3_dp, tol, 'total debt/equity')
   call assert_close(financial_leverage(16000.0_dp, 20000.0_dp), 1.25_dp, tol, 'financial leverage')
   call assert_close(gross_profit_margin(1000.0_dp, 20000.0_dp), 0.05_dp, tol, 'gpm')
   call assert_close(net_profit_margin(8000.0_dp, 20000.0_dp), 0.4_dp, tol, 'npm')
   call assert_close(ear2bey(0.08_dp), ear_to_bey(0.08_dp), tol, 'ear2bey alias')
   call assert_close(bdy2mmy(0.045_dp, 120.0_dp), bdy_to_mmy(0.045_dp, 120.0_dp), tol, 'bdy2mmy alias')
   call assert_close(gpm(1000.0_dp, 20000.0_dp), gross_profit_margin(1000.0_dp, 20000.0_dp), tol, 'gpm alias')
   call assert_close(npm(8000.0_dp, 20000.0_dp), net_profit_margin(8000.0_dp, 20000.0_dp), tol, 'npm alias')
   call assert_close(lt_d2e(8000.0_dp, 20000.0_dp), 0.4_dp, tol, 'lt_d2e alias')
   call assert_close(total_d2e(6000.0_dp, 20000.0_dp), 0.3_dp, tol, 'total_d2e alias')
   print '(a)', 'test_rates_ratios: PASS'
contains
   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: message
      if (abs(actual - expected) > tolerance) then
         print *, trim(message), actual, expected, abs(actual - expected)
         error stop 1
      end if
   end subroutine assert_close
end program test_rates_ratios
