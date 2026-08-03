program test_fcl
   use fcl
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none
   integer :: failures
   failures = 0

   call test_dates(failures)
   call test_xirr(failures)
   call test_bonds(failures)
   call test_returns(failures)

   if (failures /= 0) then
      print '(a,i0)', 'FAILED tests: ', failures
      error stop 1
   end if
   print '(a)', 'All fcl tests passed.'

contains

   subroutine check_close(actual, expected, tolerance, label, failures)
      real(dp), intent(in) :: actual, expected, tolerance
      character(*), intent(in) :: label
      integer, intent(inout) :: failures
      if (abs(actual - expected) > tolerance) then
         print '(a,2(1x,es20.10))', 'FAIL '//trim(label)//':', actual, expected
         failures = failures + 1
      end if
   end subroutine check_close

   subroutine check_true(condition, label, failures)
      logical, intent(in) :: condition
      character(*), intent(in) :: label
      integer, intent(inout) :: failures
      if (.not. condition) then
         print '(a)', 'FAIL '//trim(label)
         failures = failures + 1
      end if
   end subroutine check_true

   subroutine test_dates(failures)
      integer, intent(inout) :: failures
      type(date_type) :: d
      d = add_months(make_date(2020, 12, 31), 2)
      call check_true(d%year == 2021 .and. d%month == 2 .and. d%day == 28, &
         'add_months month end', failures)
      call check_true(days_between(make_date(2021, 1, 1), make_date(2020, 1, 1)) == 366, &
         'days_between leap year', failures)
   end subroutine test_dates

   subroutine test_xirr(failures)
      integer, intent(inout) :: failures
      real(dp) :: rate
      integer :: stat
      type(date_type) :: dates(2)
      dates = [make_date(2021, 1, 1), make_date(2022, 1, 1)]
      rate = xirr([-100.0_dp, 105.0_dp], dates, status=stat)
      call check_true(stat == 0, 'xirr status', failures)
      call check_close(rate, 0.05_dp, 1.0e-7_dp, 'xirr annual', failures)
   end subroutine test_xirr

   subroutine test_bonds(failures)
      integer, intent(inout) :: failures
      type(fixed_bond_type) :: bond
      type(bond_value_type) :: value
      type(cashflow_type) :: cf

      bond%value_date = make_date(2021, 1, 1)
      bond%maturity_date = make_date(2025, 1, 1)
      bond%redemption_value = 100.0_dp
      bond%coupon_rate = 0.05_dp
      bond%coupon_frequency = 0
      value = bond%value(make_date(2022, 1, 1), 100.0_dp)
      call check_close(value%ytm, 0.0455272763905981_dp, 2.0e-7_dp, &
         'zero coupon ytm', failures)
      call check_close(value%macaulay_duration, 3.0_dp, 1.0e-8_dp, &
         'zero coupon macaulay', failures)
      call check_close(value%modified_duration, 2.86936559941372_dp, 2.0e-6_dp, &
         'zero coupon modified', failures)

      bond%value_date = make_date(2021, 2, 1)
      bond%maturity_date = make_date(2030, 2, 1)
      bond%coupon_rate = 0.03_dp
      bond%coupon_frequency = 1
      value = bond%value(make_date(2022, 2, 1), 100.0_dp)
      call check_close(value%ytm, 0.03_dp, 2.0e-7_dp, 'coupon bond ytm', failures)
      call check_close(value%macaulay_duration, 7.23028295522156_dp, 2.0e-7_dp, &
         'coupon bond macaulay', failures)
      call check_close(value%modified_duration, 7.01969218987131_dp, 2.0e-6_dp, &
         'coupon bond modified', failures)
      cf = bond%coupon_cashflow()
      call check_true(size(cf%dates) == 9, 'coupon cashflow size', failures)
      call check_close(cf%values(1), 3.0_dp, 1.0e-12_dp, 'coupon cashflow value', failures)
   end subroutine test_bonds

   subroutine test_returns(failures)
      integer, intent(inout) :: failures
      integer :: i
      type(return_series_type) :: rtn
      real(dp), allocatable :: daily(:), cumulative(:), avc(:), dietz_values(:)

      call rtn%initialize([0, 4, 9], [100.0_dp, 103.0_dp, 110.0_dp], &
         [0.0_dp, 3.0_dp, 7.0_dp])
      call check_true(rtn%status == 0, 'return initialization', failures)
      daily = rtn%twrr_daily(1, 9)
      cumulative = rtn%twrr_cumulative(1, 9)
      avc = rtn%dietz_average_capital(1, 9)
      dietz_values = rtn%dietz(1, 9)
      call check_close(cumulative(size(cumulative)), 0.1_dp, 1.0e-12_dp, &
         'ending cumulative return', failures)
      call check_close(dietz_values(size(dietz_values)), 0.1_dp, 1.0e-12_dp, &
         'ending dietz return', failures)
      call check_close(avc(size(avc)), 100.0_dp, 1.0e-12_dp, &
         'ending average capital', failures)
      call check_true(size(daily) == 9, 'daily return length', failures)

      call rtn%initialize([1, 3, 4, 5], [100.0_dp, 102.0_dp, 103.0_dp, 104.0_dp], &
         [0.0_dp, 2.0_dp, 1.0_dp, 1.0_dp])
      cumulative = rtn%twrr_cumulative(1, 5)
      call check_true(all([(ieee_is_nan(cumulative(i)), i = 1, size(cumulative))]), &
         'first date cumulative returns are NaN', failures)
   end subroutine test_returns

end program test_fcl
