! SPDX-License-Identifier: GPL-2.0-only
! Based on FinancialMath 0.1.1, Copyright (C) 2016 Kameron Penn and Jack Schmidt.
program cashflows_and_loans
   use financialmath
   implicit none
   type(amortization_result_t) :: schedule
   real(dp) :: cash(3), times(3), rate
   logical :: ok

   cash = [300.0_dp, 400.0_dp, 500.0_dp]
   times = [1.0_dp, 2.0_dp, 3.0_dp]
   rate = irr(1000.0_dp, cash, times, ok)
   if (ok) print '(a,f10.6)', 'IRR: ', rate

   schedule = amort_table(250000.0_dp, 360.0_dp, 0.0_dp, 0.06_dp, &
      12.0_dp, 12.0_dp, 'payment')
   if (schedule%status%ok) then
      print '(a,f12.2)', 'Monthly payment: ', schedule%payment
      print '(a,f12.2)', 'Total interest:  ', schedule%total_interest
   end if
end program cashflows_and_loans
