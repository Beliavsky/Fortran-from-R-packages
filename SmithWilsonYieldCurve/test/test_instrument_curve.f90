program test_instrument_curve
   use smith_wilson_kinds, only : dp
   use smith_wilson, only : fit_smith_wilson_curve_to_instruments, make_market_instrument, &
                            market_instrument, smith_wilson_curve, sw_success
   implicit none

   type(market_instrument) :: instruments(5)
   type(smith_wilson_curve) :: curve
   real(dp), allocatable :: repriced(:), spot(:)
   real(dp) :: query(4)
   integer :: failures, info
   character(len=256) :: message

   failures = 0

   call make_market_instrument('LIBOR', 1.0_dp, 0.01_dp, instruments(1), info, message)
   call make_market_instrument('SWAP', 2.0_dp, 0.02_dp, instruments(2), info, message, frequency=1.0_dp)
   call make_market_instrument('SWAP', 3.0_dp, 0.026_dp, instruments(3), info, message, frequency=1.0_dp)
   call make_market_instrument('SWAP', 5.0_dp, 0.034_dp, instruments(4), info, message, frequency=1.0_dp)
   call make_market_instrument('BOND', 7.25_dp, 0.04_dp, instruments(5), info, message, &
                               frequency=2.0_dp, price=0.94_dp)

   call fit_smith_wilson_curve_to_instruments(instruments, 0.042_dp, 0.1_dp, curve, info, message)
   call check(info == sw_success, 'instrument fit status: '//trim(message))
   call check(curve%fitted, 'instrument curve fitted')

   repriced = curve%repriced_values()
   call check(maxval(abs(repriced - [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.94_dp])) < 2.0e-11_dp, &
              'mixed instruments reprice exactly')

   query = [1.0_dp, 10.0_dp, 50.0_dp, 200.0_dp]
   call check(all(curve%discount(query) > 0.0_dp), 'positive fitted discounts')
   spot = curve%continuous_spot(query)
   call check(abs(spot(4) - 0.042_dp) < 0.01_dp, 'long-term spot approaches UFR')
   call check_close(curve%discount(0.0_dp), 1.0_dp, 1.0e-14_dp, 'discount at zero')

   if (failures /= 0) error stop 1
   print '(a)', 'test_instrument_curve: PASS'

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

end program test_instrument_curve
