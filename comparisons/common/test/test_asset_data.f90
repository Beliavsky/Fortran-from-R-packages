program test_asset_data
   use comparison_asset_data, only : dp, asset_price_data, asset_return_data, &
      read_asset_prices, simple_returns, log_returns, asset_data_ok
   use comparison_date, only : date, date_from_iso, day_of_week, day_of_year, &
      quarter, is_month_end, operator(==)
   implicit none

   type(asset_price_data) :: prices
   type(asset_return_data) :: simple, logarithmic
   type(date) :: leap_day
   integer :: status, failures
   character(len=256) :: message
   real(dp), parameter :: tolerance = 1.0e-14_dp

   failures = 0
   call read_asset_prices('../../asset_class_etf_prices.csv', prices, status, message)
   call check(status == asset_data_ok, 'read asset prices: '//trim(message))
   if (status /= asset_data_ok) stop 1

   call check(size(prices%prices, 1) == 4607, 'price observation count')
   call check(size(prices%prices, 2) == 9, 'asset count')
   call check(prices%dates(1) == date(2007, 12, 19), 'first date')
   call check(prices%dates(4607) == date(2026, 4, 14), 'last date')
   call check(trim(prices%symbols(1)) == 'SPY', 'first symbol')
   call check(trim(prices%symbols(9)) == 'USO', 'last symbol')

   simple = simple_returns(prices)
   logarithmic = log_returns(prices)
   call check(size(simple%returns, 1) == 4606, 'simple-return observation count')
   call check(simple%dates(1) == date(2007, 12, 20), 'first return date')
   call check(logarithmic%dates(4606) == prices%dates(4607), 'last return date')
   call check(abs(simple%returns(1, 1) - &
      (104.2776_dp/103.6241_dp - 1.0_dp)) < tolerance, 'first SPY simple return')
   call check(abs(logarithmic%returns(1, 1) - &
      log(104.2776_dp/103.6241_dp)) < tolerance, 'first SPY log return')

   leap_day = date_from_iso('2024-02-29')
   call check(leap_day == date(2024, 2, 29), 'leap-day parsing')
   call check(day_of_week(date(1970, 1, 1)) == 4, 'ISO weekday')
   call check(day_of_year(leap_day) == 60, 'leap-year day number')
   call check(quarter(date(2026, 4, 14)) == 2, 'calendar quarter')
   call check(is_month_end(date(2024, 2, 29)), 'month end')

   if (failures > 0) then
      write(*, '(i0,a)') failures, ' comparison-data test(s) failed.'
      stop 1
   end if
   write(*, '(a)') 'All comparison-data tests passed.'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         failures = failures + 1
         write(*, '(a)') 'FAIL: '//trim(label)
      end if
   end subroutine check

end program test_asset_data
