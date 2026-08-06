program test_errors
   use smith_wilson_kinds, only : dp
   use smith_wilson, only : fit_smith_wilson_curve, make_market_instrument, &
                            market_instrument, smith_wilson_curve, sw_invalid_argument, &
                            sw_singular_system, sw_success, sw_unknown_instrument
   implicit none

   type(smith_wilson_curve) :: curve
   type(market_instrument) :: instrument
   real(dp) :: times(2), cashflows(2, 2), market_values(2)
   integer :: failures, info
   character(len=256) :: message

   failures = 0
   times = [1.0_dp, 2.0_dp]
   cashflows = 0.0_dp
   cashflows(1, 1) = 1.0_dp
   cashflows(2, 2) = 1.0_dp
   market_values = [0.98_dp, 0.94_dp]

   call fit_smith_wilson_curve(times, cashflows, market_values, 0.04_dp, 0.0_dp, &
                               curve, info, message)
   call check(info == sw_invalid_argument, 'reject alpha <= 0')

   cashflows(2, :) = cashflows(1, :)
   call fit_smith_wilson_curve(times, cashflows, market_values, 0.04_dp, 0.1_dp, &
                               curve, info, message)
   call check(info == sw_singular_system, 'detect singular calibration')

   call make_market_instrument('FRA', 1.0_dp, 0.02_dp, instrument, info, message)
   call check(info == sw_unknown_instrument, 'reject unknown instrument type')

   call make_market_instrument('SWAP', 1.1_dp, 0.02_dp, instrument, info, message, &
                               frequency=2.0_dp)
   call check(info == sw_invalid_argument, 'reject fractional swap payment count')

   cashflows = 0.0_dp
   cashflows(1, 1) = 1.0_dp
   cashflows(2, 2) = 1.0_dp
   call fit_smith_wilson_curve(times, cashflows, market_values, -0.01_dp, 0.1_dp, &
                               curve, info, message)
   call check(info == sw_success, 'negative UFR follows upstream warning behavior')
   call check(index(message, 'warning') > 0, 'negative UFR warning message')

   if (failures /= 0) error stop 1
   print '(a)', 'test_errors: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         failures = failures + 1
         print '(a)', 'FAIL: '//trim(label)
      end if
   end subroutine check

end program test_errors
