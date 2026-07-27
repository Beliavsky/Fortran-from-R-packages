! SPDX-License-Identifier: GPL-2.0-only
! Based on FinancialMath 0.1.1, Copyright (C) 2016 Kameron Penn and Jack Schmidt.
program test_cashflows
   use financialmath
   use test_support
   implicit none
   type(rate_conversion_t) :: rc
   type(tvm_result_t) :: tv
   type(cashflow_analysis_t) :: ca
   real(dp) :: cash(2), times(2), value, rate
   logical :: ok

   rc = rate_conv(0.12_dp, 12.0_dp, 'interest', 4.0_dp)
   call assert_close(rc%effective_interest, (1.01_dp**12)-1.0_dp)
   call assert_close(rc%effective_discount, rc%effective_interest/(1.0_dp+rc%effective_interest))
   call assert_close(rc%force, log(1.0_dp+rc%effective_interest))

   tv = solve_tvm(100.0_dp, 0.0_dp, 10.0_dp, 0.05_dp, 1.0_dp, 'fv')
   call assert_true(tv%status%ok, 'TVM FV failed')
   call assert_close(tv%future_value, 162.8894626777442_dp)
   tv = solve_tvm(100.0_dp, 162.8894626777442_dp, 10.0_dp, 0.0_dp, 1.0_dp, 'rate')
   call assert_true(tv%status%ok, 'TVM rate failed')
   call assert_close(tv%effective_rate, 0.05_dp)

   cash = [60.0_dp, 60.0_dp]
   times = [1.0_dp, 2.0_dp]
   value = npv(100.0_dp, cash, times, 0.10_dp)
   call assert_close(value, 4.13223140495867_dp)

   cash = [110.0_dp, 0.0_dp]
   times = [1.0_dp, 2.0_dp]
   rate = irr(100.0_dp, cash, times, ok)
   call assert_true(ok, 'IRR root failed')
   call assert_close(rate, 0.10_dp, 1.0e-10_dp, 1.0e-10_dp)

   cash = [5.0_dp, 105.0_dp]
   times = [1.0_dp, 2.0_dp]
   ca = cf_analysis(cash, times, 0.05_dp)
   call assert_true(ca%status%ok, 'cash-flow analysis failed')
   call assert_close(ca%present_value, 100.0_dp)
   call assert_close(ca%macaulay_duration, 1.9523809523809523_dp)
   call assert_close(ca%modified_duration, ca%macaulay_duration/1.05_dp)

   value = swap_rate([0.02_dp, 0.025_dp, 0.03_dp], 'spot_rate', ok)
   call assert_true(ok .and. value > 0.0_dp, 'swap rate failed')
   value = swap_commodity([100.0_dp, 102.0_dp], [0.02_dp, 0.03_dp], 'spot_rate', ok)
   call assert_true(ok .and. value > 100.0_dp .and. value < 102.0_dp, 'commodity swap failed')

   value = yield_dollar([10.0_dp], [0.5_dp], 100.0_dp, 120.0_dp, 1.0_dp)
   call assert_close(value, 10.0_dp/105.0_dp)
   value = yield_time([5.0_dp], [100.0_dp, 110.0_dp])
   call assert_close(value, 110.0_dp/105.0_dp-1.0_dp)

   print '(a)', 'test_cashflows: PASS'
end program test_cashflows
