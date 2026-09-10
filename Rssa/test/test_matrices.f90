! SPDX-License-Identifier: GPL-2.0-or-later
! Test/example code for the modern Fortran translation of the upstream Rssa package.
program test_matrices
   use rssa, only : dp, hankel_matrix, hankelize_matrix, lag_covariance, trajectory_2d
   implicit none
   real(dp) :: x(5), field(3, 4)
   real(dp), allocatable :: h(:, :), y(:), t(:, :), lag(:)
   integer :: i
   x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   h = hankel_matrix(x, 3)
   call assert_true(all(shape(h) == [3, 3]), 'Hankel shape')
   call assert_close(maxval(abs(h(:, 1) - [1.0_dp, 2.0_dp, 3.0_dp])), 0.0_dp, 1.0e-12_dp, 'Hankel values')
   y = hankelize_matrix(h)
   call assert_close(maxval(abs(y - x)), 0.0_dp, 1.0e-12_dp, 'Hankelization')
   field = reshape([(real(i, dp), i=1, 12)], [3, 4])
   t = trajectory_2d(field, [2, 2])
   call assert_true(all(shape(t) == [4, 6]), '2-D trajectory shape')
   lag = lag_covariance(x, 3, .false.)
   call assert_true(size(lag) == 3, 'lag covariance size')
   print *, 'test_matrices: PASS'
contains
   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Condition required to be true.
      character(len=*), intent(in) :: message !! Diagnostic label printed on failure.
      if (.not. condition) error stop message
   end subroutine assert_true
   subroutine assert_close(value, expected, tolerance, message)
      real(dp), intent(in) :: value !! Computed scalar value.
      real(dp), intent(in) :: expected !! Expected scalar reference.
      real(dp), intent(in) :: tolerance !! Maximum accepted absolute error.
      character(len=*), intent(in) :: message !! Diagnostic label printed on failure.
      if (abs(value - expected) > tolerance) error stop message
   end subroutine assert_close
end program test_matrices
