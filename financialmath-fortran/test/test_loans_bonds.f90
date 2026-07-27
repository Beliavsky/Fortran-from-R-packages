! SPDX-License-Identifier: GPL-2.0-only
! Based on FinancialMath 0.1.1, Copyright (C) 2016 Kameron Penn and Jack Schmidt.
program test_loans_bonds
   use financialmath
   use test_support
   implicit none
   type(amortization_result_t) :: am
   type(amort_period_result_t) :: ap
   type(bond_result_t) :: br
   real(dp) :: expected_payment

   expected_payment = 100000.0_dp*(0.06_dp/12.0_dp)/(1.0_dp-(1.0_dp+0.06_dp/12.0_dp)**(-360.0_dp))
   am = amort_table(100000.0_dp, 360.0_dp, 0.0_dp, 0.06_dp, 12.0_dp, 12.0_dp, 'payment')
   call assert_true(am%status%ok, 'amortization table failed')
   call assert_close(am%payment, expected_payment)
   call assert_close(am%balance(size(am%balance)), 0.0_dp, 1.0e-7_dp, 0.0_dp)
   call assert_close(sum(am%principal), 100000.0_dp, 1.0e-6_dp, 1.0e-10_dp)

   ap = amort_period(100000.0_dp, 360.0_dp, 0.0_dp, 0.06_dp, 12.0_dp, 12.0_dp, 1.0_dp, 'payment')
   call assert_true(ap%status%ok, 'amortization period failed')
   call assert_close(ap%payment, expected_payment)
   call assert_close(ap%interest_paid, 500.0_dp)
   call assert_close(ap%principal_paid, expected_payment-500.0_dp)

   br = bond(1000.0_dp, 0.05_dp, 1000.0_dp, 10.0_dp, 0.04_dp, 1.0_dp, 1.0_dp, 0.0_dp, .false.)
   call assert_true(br%status%ok, 'bond calculation failed')
   call assert_close(br%price, 1081.108957793552_dp)
   call assert_true(br%premium > 0.0_dp .and. br%discount < 1.0e-12_dp, 'bond premium failed')
   call assert_true(br%macaulay_duration > br%modified_duration, 'bond duration ordering failed')

   print '(a)', 'test_loans_bonds: PASS'
end program test_loans_bonds
