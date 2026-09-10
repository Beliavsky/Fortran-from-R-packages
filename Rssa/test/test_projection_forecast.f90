! SPDX-License-Identifier: GPL-2.0-or-later
! Test/example code for the modern Fortran translation of the upstream Rssa package.
program test_projection_forecast
   use rssa_forecast, only : rforecast_pssa, vforecast_pssa
   use rssa_kinds, only : dp
   use rssa_projection, only : decompose_pssa
   use rssa_types, only : rssa_success, ssa_result
   implicit none
   type(ssa_result) :: object, projected
   real(dp), allocatable :: recurrent(:), vector(:)
   real(dp) :: column_projector(3, 1), row_projector(4, 1)
   real(dp) :: series(6)
   integer :: info, n_special

   allocate(object%series(6), object%sigma(2), object%u(3, 2), object%v(4, 2))
   object%series = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
   object%window = 3
   object%n_special_right = 1
   object%n_special_left = 0
   object%sigma = [1.3_dp, 0.8_dp]
   object%u(:, 1) = [0.7_dp, 0.1_dp, 0.2_dp]
   object%u(:, 2) = [0.2_dp, 0.5_dp, 0.4_dp]
   object%v(:, 1) = [0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp]
   object%v(:, 2) = [0.1_dp, 0.4_dp, 0.2_dp, 0.5_dp]

   recurrent = rforecast_pssa(object, [1, 2], len=2, only_new=.true., info=info)
   if (info /= rssa_success) error stop "rforecast_pssa status"
   call assert_close(recurrent, [0.1620952380952381_dp, 0.1374036281179138_dp], 2.0e-12_dp, &
                     "rforecast_pssa reference")

   vector = vforecast_pssa(object, [1, 2], len=2, only_new=.true., info=info)
   if (info /= rssa_success) error stop "vforecast_pssa status"
   call assert_close(vector, [0.2661103366271679_dp, 0.2102865551315530_dp], 2.0e-11_dp, &
                     "vforecast_pssa reference")

   series = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
   column_projector(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
   row_projector(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   projected = decompose_pssa(series, 3, column_projector, row_projector, neig=0, n_special=n_special, info=info)
   if (info /= rssa_success) error stop "decompose_pssa status"
   if (n_special /= 2) error stop "projection special count"
   if (projected%n_special_right /= 1) error stop "right special metadata"
   if (projected%n_special_left /= 1) error stop "left special metadata"

   print *, "test_projection_forecast: PASS"
contains
   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:) !! Values produced by the Fortran implementation.
      real(dp), intent(in) :: expected(:) !! Independently computed reference values.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute error.
      character(len=*), intent(in) :: label !! Diagnostic label emitted on failure.

      if (size(actual) /= size(expected)) error stop label // ": size mismatch"
      if (size(actual) > 0) then
         if (maxval(abs(actual - expected)) > tolerance) error stop label // ": value mismatch"
      end if
   end subroutine assert_close
end program test_projection_forecast
