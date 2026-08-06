program test_instruments
   use smith_wilson_kinds, only : dp
   use smith_wilson, only : create_cashflow_matrix, create_market_value_vector, &
                            create_time_vector, get_instrument_cashflows, &
                            make_market_instrument, market_instrument, sw_success
   implicit none

   type(market_instrument) :: instruments(3)
   real(dp), allocatable :: times(:), schedule_times(:), cashflows(:), matrix(:, :), values(:)
   integer :: failures, info
   character(len=256) :: message

   failures = 0

   call make_market_instrument('libor', 0.5_dp, 0.04_dp, instruments(1), info, message)
   call check(info == sw_success, 'construct LIBOR')
   call make_market_instrument('SWAP', 2.0_dp, 0.03_dp, instruments(2), info, message, &
                               frequency=2.0_dp)
   call check(info == sw_success, 'construct swap')
   call make_market_instrument('BOND', 2.25_dp, 0.04_dp, instruments(3), info, message, &
                               frequency=2.0_dp, price=0.96_dp)
   call check(info == sw_success, 'construct bond')

   call get_instrument_cashflows(instruments(1), schedule_times, cashflows, info, message)
   call check(info == sw_success, 'LIBOR schedule status')
   call check(size(schedule_times) == 1, 'LIBOR schedule length')
   call check_close(schedule_times(1), 0.5_dp, 1.0e-15_dp, 'LIBOR time')
   call check_close(cashflows(1), 1.02_dp, 1.0e-15_dp, 'LIBOR cashflow')

   call get_instrument_cashflows(instruments(2), schedule_times, cashflows, info, message)
   call check(size(schedule_times) == 4, 'swap schedule length')
   call check(maxval(abs(schedule_times - [0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp])) < 1.0e-15_dp, &
              'swap schedule times')
   call check(maxval(abs(cashflows - [0.015_dp, 0.015_dp, 0.015_dp, 1.015_dp])) < 1.0e-15_dp, &
              'swap cashflows')

   call get_instrument_cashflows(instruments(3), schedule_times, cashflows, info, message)
   call check(size(schedule_times) == 5, 'bond schedule length')
   call check(maxval(abs(schedule_times - [0.25_dp, 0.75_dp, 1.25_dp, 1.75_dp, 2.25_dp])) < 1.0e-15_dp, &
              'bond schedule times')
   call check(maxval(abs(cashflows - [0.02_dp, 0.02_dp, 0.02_dp, 0.02_dp, 1.02_dp])) < 1.0e-15_dp, &
              'bond cashflows')

   call create_time_vector(instruments, times, info, message)
   call check(info == sw_success, 'combined time-vector status')
   call check(maxval(abs(times - [0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp, 1.25_dp, &
                                  1.5_dp, 1.75_dp, 2.0_dp, 2.25_dp])) < 1.0e-15_dp, &
              'combined unique times')

   call create_cashflow_matrix(instruments, times, matrix, info, message)
   call check(info == sw_success, 'cashflow-matrix status')
   call check(size(matrix, 1) == 3 .and. size(matrix, 2) == 9, 'cashflow-matrix shape')
   call check_close(matrix(1, 2), 1.02_dp, 1.0e-15_dp, 'LIBOR placement')
   call check_close(matrix(2, 8), 1.015_dp, 1.0e-15_dp, 'swap final placement')
   call check_close(matrix(3, 9), 1.02_dp, 1.0e-15_dp, 'bond final placement')

   call create_market_value_vector(instruments, values, info, message)
   call check(maxval(abs(values - [1.0_dp, 1.0_dp, 0.96_dp])) < 1.0e-15_dp, &
              'market values')

   if (failures /= 0) error stop 1
   print '(a)', 'test_instruments: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         failures = failures + 1
         print '(a)', 'FAIL: '//trim(label)
      end if
   end subroutine check

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label

      call check(abs(actual - expected) <= tolerance, label)
   end subroutine check_close

end program test_instruments
