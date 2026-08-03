! SPDX-License-Identifier: GPL-3.0-or-later
program test_swaps
   use quant_bond_curves
   implicit none
   type(qbc_swap) :: swap
   type(qbc_curve) :: curve
   real(dp) :: value
   integer :: st

   allocate(curve%terms(3), curve%rates(3))
   curve%terms = [0.0_dp, 2.0_dp, 5.0_dp]
   curve%rates = 0.04_dp
   curve%approximation = 2
   curve%rate_type = qbc_rate_continuous
   swap%analysis_date = make_date(2024,1,1)
   swap%maturity = make_date(2029,1,1)
   swap%frequency = 1
   swap%rate_type = qbc_rate_continuous
   swap%coupon_rate1 = 0.05_dp
   swap%coupon_rate2 = 0.03_dp
   swap%principal1 = 100.0_dp
   swap%principal2 = 100.0_dp
   swap%exchange_rate = 1.0_dp
   swap%legs = qbc_leg_fixed_fixed
   value = valuation_swaps(swap, curve, curve, status=st)
   call assert_true(st == qbc_success .and. value > 0.0_dp, 'higher fixed coupon has positive value')
   swap%coupon_rate2 = swap%coupon_rate1
   call assert_close(valuation_swaps(swap, curve, curve), 0.0_dp, 1.0e-12_dp, 'identical legs cancel')

   print '(a)', 'test_swaps: PASS'
contains
   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) error stop label
   end subroutine
   subroutine assert_close(x, y, tol, label)
      real(dp), intent(in) :: x, y, tol
      character(len=*), intent(in) :: label
      if (abs(x-y) > tol) then
         write(*,*) label, x, y
         error stop 'assert_close'
      end if
   end subroutine
end program test_swaps
