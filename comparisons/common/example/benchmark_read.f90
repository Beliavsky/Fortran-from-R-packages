program benchmark_read
   use iso_fortran_env, only : int64
   use comparison_asset_data, only : dp, asset_price_data, asset_return_data, &
      read_asset_prices, read_asset_prices_binary, log_returns, asset_data_ok
   implicit none

   integer, parameter :: repetitions = 20
   type(asset_price_data) :: prices
   type(asset_return_data) :: returns
   integer :: status, iteration
   integer(int64) :: clock_rate, start_clock, end_clock
   real(dp) :: first_seconds, repeated_seconds, binary_first_seconds
   real(dp) :: binary_repeated_seconds, returns_seconds, checksum
   character(len=512) :: filename, binary_filename
   character(len=256) :: message

   call get_command_argument(1, filename)
   call get_command_argument(2, binary_filename)
   if (len_trim(filename) == 0) filename = '../../asset_class_etf_prices.csv'
   if (len_trim(binary_filename) == 0) binary_filename = '../../asset_class_etf_prices.bin'
   call system_clock(count_rate=clock_rate)

   call system_clock(start_clock)
   call read_asset_prices(trim(filename), prices, status, message)
   call system_clock(end_clock)
   if (status /= asset_data_ok) error stop trim(message)
   first_seconds = elapsed_seconds(start_clock, end_clock, clock_rate)

   call system_clock(start_clock)
   do iteration = 1, repetitions
      call read_asset_prices(trim(filename), prices, status, message)
      if (status /= asset_data_ok) error stop trim(message)
   end do
   call system_clock(end_clock)
   repeated_seconds = elapsed_seconds(start_clock, end_clock, clock_rate)

   call system_clock(start_clock)
   call read_asset_prices_binary(trim(binary_filename), prices, status, message)
   call system_clock(end_clock)
   if (status /= asset_data_ok) error stop trim(message)
   binary_first_seconds = elapsed_seconds(start_clock, end_clock, clock_rate)

   call system_clock(start_clock)
   do iteration = 1, repetitions
      call read_asset_prices_binary(trim(binary_filename), prices, status, message)
      if (status /= asset_data_ok) error stop trim(message)
   end do
   call system_clock(end_clock)
   binary_repeated_seconds = elapsed_seconds(start_clock, end_clock, clock_rate)

   checksum = 0.0_dp
   call system_clock(start_clock)
   do iteration = 1, repetitions
      returns = log_returns(prices)
      checksum = checksum + sum(returns%returns)
   end do
   call system_clock(end_clock)
   returns_seconds = elapsed_seconds(start_clock, end_clock, clock_rate)

   write(*, '(a,i0)') 'Price rows: ', size(prices%prices, 1)
   write(*, '(a,i0)') 'Assets: ', size(prices%prices, 2)
   write(*, '(a,f10.6,a)') 'First CSV read: ', first_seconds, ' seconds'
   write(*, '(i0,a,f10.6,a)') repetitions, ' cached CSV reads: ', &
      repeated_seconds/real(repetitions, dp), ' seconds/read'
   write(*, '(a,f10.6,a)') 'First binary read: ', binary_first_seconds, ' seconds'
   write(*, '(i0,a,f10.6,a)') repetitions, ' cached binary reads: ', &
      binary_repeated_seconds/real(repetitions, dp), ' seconds/read'
   write(*, '(i0,a,f10.6,a)') repetitions, ' return constructions: ', &
      returns_seconds/real(repetitions, dp), ' seconds/construction'
   write(*, '(a,es16.8)') 'Checksum: ', checksum

contains

   pure real(dp) function elapsed_seconds(first, last, rate)
      integer(int64), intent(in) :: first, last, rate

      elapsed_seconds = real(last - first, dp)/real(rate, dp)
   end function elapsed_seconds

end program benchmark_read
