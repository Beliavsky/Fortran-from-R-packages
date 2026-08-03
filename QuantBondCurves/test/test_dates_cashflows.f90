! SPDX-License-Identifier: GPL-3.0-or-later
program test_dates_cashflows
   use quant_bond_curves
   implicit none
   type(qbc_coupon_schedule) :: schedule
   type(qbc_date) :: d
   real(dp), allocatable :: cf(:)
   integer :: st

   call parse_date('2024-02-29', d, st)
   call assert_true(st == qbc_success .and. date_string(d) == '2024-02-29', 'date parsing')
   call assert_close(discount_time(make_date(2024,2,29), make_date(2025,2,28)), 1.0_dp, 1.0e-14_dp, 'leap anniversary')
   call assert_close(discount_time(make_date(2024,2,29), make_date(2028,2,29)), 4.0_dp, 1.0e-14_dp, 'four-year anniversary')
   call assert_close(year_fraction(make_date(2024,1,1), make_date(2024,7,1), 'ACT/366'), &
                     discount_time(make_date(2024,1,1), make_date(2024,7,1)), 1.0e-14_dp, 'fallback convention')
   call assert_close(year_fraction(make_date(2024,1,1), make_date(2024,7,1), 'ACT/365L'), 182.0_dp/366.0_dp, 1.0e-14_dp, 'ACT/365L')

   call coupon_dates(make_date(2026,5,29), make_date(2020,2,29), 2, qbc_schedule_short_last, &
                     'F', schedule, trade_date=make_date(2020,2,29), status=st)
   call assert_true(st == qbc_success .and. size(schedule%dates) == 13, 'short-last schedule size')
   call assert_true(date_string(schedule%dates(1)) == '2020-08-29', 'first short-last coupon')
   call assert_true(date_string(schedule%dates(13)) == '2026-05-29', 'last short-last coupon')
   call assert_true(date_string(schedule%effective_dates(1)) == '2020-08-31', 'following convention')

   call coupon_cashflows([make_date(2024,7,1), make_date(2025,1,1)], [0.06_dp], 100.0_dp, &
                         'ACT/365', cf, make_date(2024,1,1), .false., st)
   call assert_true(st == qbc_success .and. size(cf) == 2, 'cash-flow size')
   call assert_close(cf(1), 0.06_dp*100.0_dp*182.0_dp/365.0_dp, 1.0e-12_dp, 'first coupon')
   call assert_close(cf(2), 100.0_dp + 0.06_dp*100.0_dp*184.0_dp/365.0_dp, 1.0e-12_dp, 'redemption coupon')

   print '(a)', 'test_dates_cashflows: PASS'
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
end program test_dates_cashflows
