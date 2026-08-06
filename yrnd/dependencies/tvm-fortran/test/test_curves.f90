! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Juan Manuel Truppia
program test_curves
   use tvm
   implicit none
   real(dp), parameter :: tol = 2.0e-8_dp
   real(dp) :: flat(12), discount(12), upward(3), expected(3), swap_rates(3)

   call one_rate_tests()

   flat = 0.1_dp
   discount = fut_to_disc(flat)
   call assert_array_close(discount, 1.0_dp / cumulative_product(1.0_dp + flat), tol, "flat discount")
   call assert_array_close(disc_to_zero_eff(discount), flat, tol, "flat zero")
   call assert_array_close(disc_to_french(discount), flat, tol, "flat french")
   call assert_array_close(disc_to_german(discount), flat, tol, "flat german")
   call assert_array_close(disc_to_fut(discount), flat, tol, "flat fut")
   call assert_array_close(disc_to_swap(discount), flat, tol, "flat swap")

   upward = [0.1_dp, 0.2_dp, 0.3_dp]
   discount(1:3) = fut_to_disc(upward)
   expected = cumulative_product(1.0_dp + upward) ** [1.0_dp, 0.5_dp, 1.0_dp / 3.0_dp] - 1.0_dp
   call assert_array_close(disc_to_fut(discount(1:3)), upward, tol, "upward fut")
   call assert_array_close(disc_to_zero_eff(discount(1:3)), expected, tol, "upward zero")
   swap_rates = disc_to_swap(discount(1:3))
   call assert_close(find_rate(3, discount(1:3), "bullet"), swap_rates(3), tol, "bullet rate")

   call round_trip_tests()
   print '(a)', "test_curves: PASS"

contains

   subroutine one_rate_tests()
      real(dp) :: r(1), d(1)
      r = 0.1_dp
      d = zero_eff_to_disc(r)
      call assert_close(d(1), 1.0_dp / 1.1_dp, tol, "one discount")
      call assert_array_close(disc_to_zero_eff(d), r, tol, "one zero")
      call assert_array_close(disc_to_french(d), r, tol, "one french")
      call assert_array_close(disc_to_german(d), r, tol, "one german")
      call assert_array_close(disc_to_fut(d), r, tol, "one fut")
      call assert_array_close(disc_to_swap(d), r, tol, "one swap")
   end subroutine one_rate_tests

   subroutine round_trip_tests()
      real(dp) :: z(4), d(4)
      z = [0.01_dp, 0.015_dp, 0.02_dp, 0.025_dp]
      d = zero_nom_to_disc(z)
      call assert_array_close(disc_to_zero_nom(d), z, tol, "nominal round trip")
      d = zero_cont_to_disc(z)
      call assert_array_close(disc_to_zero_cont(d), z, tol, "continuous round trip")
      call assert_array_close(dir_to_eff(eff_to_dir(z)), z, tol, "direct round trip")
   end subroutine round_trip_tests

   function cumulative_product(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x))
      integer :: i
      if (size(x) == 0) return
      y(1) = x(1)
      do i = 2, size(x)
         y(i) = y(i - 1) * x(i)
      end do
   end function cumulative_product

   subroutine assert_close(actual, expected_value, tolerance, label)
      real(dp), intent(in) :: actual, expected_value, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual - expected_value) > tolerance) then
         print '(a,2(1x,es24.16))', trim(label), actual, expected_value
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_array_close(actual, expected_value, tolerance, label)
      real(dp), intent(in) :: actual(:), expected_value(:), tolerance
      character(len=*), intent(in) :: label
      if (size(actual) /= size(expected_value) .or. any(abs(actual - expected_value) > tolerance)) then
         print '(a)', trim(label)
         error stop 1
      end if
   end subroutine assert_array_close

end program test_curves
