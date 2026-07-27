! SPDX-License-Identifier: GPL-2.0-only
! Based on FinancialMath 0.1.1, Copyright (C) 2016 Kameron Penn and Jack Schmidt.
program options_and_forwards
   use financialmath
   implicit none
   type(option_order1_t) :: bs
   type(payoff_table_t) :: spread
   type(forward_result_t) :: fw

   bs = bls_order1(100.0_dp, 100.0_dp, 0.04_dp, 0.5_dp, 0.25_dp, 0.01_dp)
   spread = bull_call_bls(100.0_dp, 90.0_dp, 110.0_dp, 0.04_dp, 0.5_dp, 0.25_dp)
   fw = forward_contract(100.0_dp, 0.5_dp, 0.04_dp, 'long', &
      'continuous', 0.0_dp, 1.0_dp, 0.01_dp, -1.0_dp)

   print '(a,f10.4)', 'Call price:       ', bs%call_price
   print '(a,f10.4)', 'Call delta:       ', bs%call_delta
   print '(a,f10.4)', 'Bull spread cost: ', spread%premiums(3)
   print '(a,f10.4)', 'Forward price:    ', fw%delivery_price
end program options_and_forwards
