! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Juan Manuel Truppia
program test_rate_curve
   use tvm
   implicit none
   type(rate_curve_t) :: curve, callback_curve
   real(dp) :: x(3), values(3), pv

   curve = rate_curve_from_rates([0.1_dp, 0.2_dp, 0.3_dp], "zero_eff")
   call assert_array_close(curve%rate_grid("zero_eff"), [0.1_dp, 0.2_dp, 0.3_dp], 2.0e-8_dp, "curve rates")
   x = [0.0_dp, 1.0_dp, 2.0_dp]
   values = curve%discount(x)
   call assert_close(values(1), 1.0_dp, 1.0e-14_dp, "discount zero")
   call assert_close(values(2), 1.0_dp / 1.1_dp, 1.0e-12_dp, "discount one")
   pv = disc_value(curve, [-1.0_dp, 1.1_dp], [0.0_dp, 1.0_dp])
   call assert_close(pv, 0.0_dp, 2.0e-12_dp, "disc value")

   callback_curve = rate_curve_from_discount_function(discount_function, [1.0_dp, 2.0_dp, 3.0_dp])
   call assert_close(callback_curve%discount(1.5_dp), 0.4_dp, 1.0e-14_dp, "callback")

   curve = rate_curve_from_rates([0.12_dp, 0.12_dp, 0.12_dp], "zero_eff", rate_scale=12.0_dp)
   call assert_array_close(curve%rate_grid("zero_eff"), [0.12_dp, 0.12_dp, 0.12_dp], 2.0e-8_dp, "scaling")

   print '(a)', "test_rate_curve: PASS"

contains

   function discount_function(time) result(value)
      real(dp), intent(in) :: time
      real(dp) :: value
      value = 1.0_dp / (1.0_dp + time)
   end function discount_function

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

end program test_rate_curve
