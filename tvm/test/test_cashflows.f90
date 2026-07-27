! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Juan Manuel Truppia
program test_cashflows
   use tvm, only : dp, adjust_disc, cft, npv, xnpv, irr, xirr, pmt, rate, loan, cashflow, rem, loan_t, &
      continuous_compounding
   implicit none
   type(loan_t) :: l
   real(dp), allocatable :: cf(:), balance(:)
   real(dp) :: value
   integer :: status

   call assert_close(npv(0.01_dp, [-1.0_dp, 0.5_dp, 0.9_dp], [0.0_dp, 1.0_dp, 3.0_dp]), &
      0.3685806380853751_dp, 1.0e-13_dp, "npv")
   call assert_close(xnpv(0.01_dp, [-1.0_dp, 0.5_dp, 0.9_dp], [0, 45, 99]), &
      0.3969613031182203_dp, 1.0e-13_dp, "xnpv dates")
   call assert_close(xnpv(0.01_dp, [-1.0_dp, 0.5_dp, 0.9_dp], &
      [0.0_dp, 45.0_dp / 365.0_dp, 99.0_dp / 365.0_dp]), &
      0.3969613031182203_dp, 1.0e-13_dp, "xnpv tau")
   call assert_close(xnpv(0.01_dp, [-1.0_dp, 1.0_dp], [0.0_dp, 1.0_dp], 0.0_dp), &
      -1.0_dp + 1.0_dp / 1.01_dp, 1.0e-14_dp, "xnpv simple")
   call assert_close(xnpv(0.01_dp, [-1.0_dp, 1.0_dp], [0.0_dp, 1.0_dp], continuous_compounding), &
      -1.0_dp + exp(-0.01_dp), 1.0e-14_dp, "xnpv continuous")

   value = irr([-1.0_dp, 0.5_dp, 0.9_dp], [0.0_dp, 1.0_dp, 3.0_dp], status=status)
   call assert_true(status == 0, "irr status")
   call assert_close(value, 0.1641203206559595_dp, 2.0e-10_dp, "irr")
   value = xirr([-1.0_dp, 1.5_dp], [0, 365], status=status)
   call assert_true(status == 0, "xirr status")
   call assert_close(value, 0.5_dp, 2.0e-10_dp, "xirr")

   call assert_close(pmt(100.0_dp, 10, 0.05_dp), 12.95045749654566_dp, 1.0e-12_dp, "pmt")
   value = rate(100.0_dp, 10, 15.0_dp, tol=1.0e-10_dp, status=status)
   call assert_true(status == 0, "rate status")
   call assert_close(pmt(100.0_dp, 10, value), 15.0_dp, 1.0e-7_dp, "rate inversion")
   call assert_close(cft(100.0_dp, 10, 0.05_dp, 1.0_dp, 0.1_dp, status), &
      0.05363357185856432_dp, 2.0e-4_dp, "cft")

   l = loan(0.1_dp, 4, 1.0_dp, "bullet")
   cf = cashflow(l)
   call assert_array_close(cf, [0.1_dp, 0.1_dp, 0.1_dp, 1.1_dp], 1.0e-14_dp, "bullet")
   l = loan(0.1_dp, 4, 1.0_dp, "german")
   call assert_array_close(l%cf, [0.35_dp, 0.325_dp, 0.3_dp, 0.275_dp], 1.0e-14_dp, "german")
   l = loan(0.1_dp, 4, 1.0_dp, "french")
   call assert_array_close(l%cf, spread(pmt(1.0_dp, 4, 0.1_dp), 1, 4), 1.0e-14_dp, "french")
   l = loan(0.1_dp, 5, 1.0_dp, "bullet", grace_int=1, grace_amort=2)
   call assert_array_close(l%cf, [0.0_dp, 0.11_dp, 0.11_dp, 0.11_dp, 1.21_dp], 1.0e-14_dp, "grace")

   balance = rem([0.4_dp, 0.4_dp, 0.4_dp, 0.4_dp], 1.0_dp, 0.2_dp)
   call assert_array_close(balance, [0.8_dp, 0.56_dp, 0.272_dp, -0.0736_dp], 1.0e-13_dp, "rem")
   call assert_array_close(adjust_disc([0.99_dp, 0.98_dp], 0.01_dp), &
      [0.9802950787206653_dp, 0.9608813534412108_dp], 1.0e-13_dp, "adjust_disc")

   print '(a)', "test_cashflows: PASS"

contains

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual - expected) > tolerance) then
         print '(a,2(1x,es24.16))', trim(label), actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_array_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      character(len=*), intent(in) :: label
      if (size(actual) /= size(expected) .or. any(abs(actual - expected) > tolerance)) then
         print '(a)', trim(label)
         error stop 1
      end if
   end subroutine assert_array_close

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         print '(a)', trim(label)
         error stop 1
      end if
   end subroutine assert_true

end program test_cashflows
