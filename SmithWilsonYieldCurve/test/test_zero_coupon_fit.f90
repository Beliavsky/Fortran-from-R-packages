program test_zero_coupon_fit
   use smith_wilson_kinds, only : dp
   use smith_wilson, only : fit_smith_wilson_curve, smith_wilson_curve, sw_success
   implicit none

   type(smith_wilson_curve) :: curve
   real(dp) :: times(2), cashflows(2, 2), market_values(2), query(6)
   real(dp), allocatable :: prices(:), repriced(:), kernel(:, :)
   integer :: failures, info
   character(len=256) :: message

   failures = 0
   times = [5.0_dp, 20.0_dp]
   cashflows = 0.0_dp
   cashflows(1, 1) = 1.0_dp
   cashflows(2, 2) = 1.0_dp
   market_values = [0.88_dp, 0.37_dp]

   call fit_smith_wilson_curve(times, cashflows, market_values, 0.042_dp, 0.1_dp, &
                               curve, info, message)
   call check(info == sw_success, 'fit status: '//trim(message))
   call check(curve%fitted, 'curve fitted flag')
   call check_close(curve%xi(1), 2.5248944821326880_dp, 2.0e-12_dp, 'first weight')
   call check_close(curve%xi(2), -1.5685333215051975_dp, 2.0e-12_dp, 'second weight')

   query = [1.0_dp, 5.0_dp, 10.0_dp, 20.0_dp, 50.0_dp, 100.0_dp]
   prices = curve%discount(query)
   call check_close(prices(1), 0.97976001432652993_dp, 2.0e-13_dp, 'P(1)')
   call check_close(prices(2), 0.88_dp, 2.0e-13_dp, 'P(5)')
   call check_close(prices(3), 0.69747140295675303_dp, 2.0e-13_dp, 'P(10)')
   call check_close(prices(4), 0.37_dp, 2.0e-13_dp, 'P(20)')
   call check_close(prices(5), 0.083071643243299936_dp, 2.0e-13_dp, 'P(50)')
   call check_close(prices(6), 0.010033214193183693_dp, 2.0e-13_dp, 'P(100)')
   call check_close(curve%discount(0.0_dp), 1.0_dp, 1.0e-15_dp, 'P(0)')

   repriced = curve%repriced_values()
   call check(maxval(abs(repriced - market_values)) < 5.0e-13_dp, 'instrument repricing')

   kernel = curve%compound_kernel(query)
   call check(size(kernel, 1) == 2 .and. size(kernel, 2) == size(query), &
              'compound kernel shape')

   if (failures /= 0) error stop 1
   print '(a)', 'test_zero_coupon_fit: PASS'

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

end program test_zero_coupon_fit
