! SPDX-License-Identifier: GPL-3.0-only
! Derived from R package ftsa 6.7 by Rob Hyndman and Han Lin Shang.
program test_statistics
   use ftsa, only : dp, ftsa_centre_mean, ftsa_difference, ftsa_error_metric, ftsa_extract_time
   use ftsa, only : ftsa_mean_coordinate, ftsa_median_coordinate, ftsa_quantiles, ftsa_variance_coordinate
   implicit none
   real(dp) :: y(2, 4), mean_curve(2), median_curve(2), var_curve(2), value
   real(dp), allocatable :: diff(:, :), extracted(:, :), q(:, :)
   integer :: status

   y = reshape([1.0_dp, 10.0_dp, 2.0_dp, 20.0_dp, 3.0_dp, 30.0_dp, 4.0_dp, 40.0_dp], [2, 4])
   call ftsa_mean_coordinate(y, mean_curve)
   call check(maxval(abs(mean_curve - [2.5_dp, 25.0_dp])) < 1.0e-12_dp, 'mean')
   call ftsa_median_coordinate(y, median_curve)
   call check(maxval(abs(median_curve - [2.5_dp, 25.0_dp])) < 1.0e-12_dp, 'median')
   call ftsa_variance_coordinate(y, var_curve)
   call check(maxval(abs(var_curve - [1.25_dp, 125.0_dp])) < 1.0e-12_dp, 'variance')
   call ftsa_centre_mean(y, mean_curve)
   call check(abs(mean_curve(1) - 2.5_dp) < 1.0e-12_dp, 'centre mean')
   call ftsa_difference(y, differenced=diff, status=status)
   call check(status == 0 .and. all(shape(diff) == [2, 3]), 'difference shape')
   call check(maxval(abs(diff(:, 1) - [1.0_dp, 10.0_dp])) < 1.0e-12_dp, 'difference values')
   call ftsa_extract_time(y, [2, 4], extracted, status)
   call check(status == 0 .and. all(shape(extracted) == [2, 2]), 'extract')
   call ftsa_quantiles(y, [0.25_dp, 0.5_dp, 0.75_dp], q, status=status)
   call check(status == 0 .and. abs(q(1, 2) - 2.5_dp) < 1.0e-12_dp, 'quantiles')
   call ftsa_error_metric([1.0_dp, 2.0_dp], [2.0_dp, 4.0_dp], 'mae', value, status=status)
   call check(status == 0 .and. abs(value - 1.5_dp) < 1.0e-12_dp, 'mae')
   call ftsa_error_metric([1.0_dp, 2.0_dp], [2.0_dp, 4.0_dp], 'rmse', value, status=status)
   call check(status == 0 .and. abs(value - 1.581139_dp) < 1.0e-6_dp, 'rmse')

contains
   subroutine check(condition, message)
      logical, intent(in) :: condition !! Assertion condition expected to be true.
      character(len=*), intent(in) :: message !! Diagnostic emitted if the assertion fails.
      if (.not. condition) error stop message
   end subroutine check
end program test_statistics
