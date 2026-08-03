! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_kernels_hac
   use sandwich, only : dp, SANDWICH_SUCCESS, kernel_weight, pava_fitted, meat_hac, &
      prewhite_var, bandwidth_andrews, bandwidth_newey_west
   implicit none

   real(dp) :: scores(5, 2), weights(2), xiso(4), series(40, 2), bw
   real(dp), allocatable :: fitted(:), m(:, :), residuals(:, :), recolor(:, :)
   integer :: status, i

   call assert_close(kernel_weight(0.5_dp, 'Bartlett'), 0.5_dp, 1.0e-14_dp, 'Bartlett')
   call assert_close(kernel_weight(0.25_dp, 'Parzen'), 0.71875_dp, 1.0e-14_dp, 'Parzen')
   call assert_close(kernel_weight(0.5_dp, 'Tukey-Hanning'), 0.5_dp, 1.0e-14_dp, 'Tukey')
   call assert_close(kernel_weight(0.0_dp, 'QS'), 1.0_dp, 1.0e-14_dp, 'QS origin')

   xiso = [3.0_dp, 1.0_dp, 2.0_dp, 4.0_dp]
   call pava_fitted(xiso, fitted, status)
   call assert_true(status == SANDWICH_SUCCESS, 'PAVA status')
   call assert_vector_close(fitted, [2.0_dp, 2.0_dp, 2.0_dp, 4.0_dp], 1.0e-14_dp, 'PAVA')

   scores = reshape([ &
      1.0_dp, -1.0_dp, 2.0_dp, 0.0_dp, 1.0_dp, &
      2.0_dp,  1.0_dp, 0.0_dp,-2.0_dp, 1.0_dp], [5, 2])
   weights = [1.0_dp, 0.5_dp]
   call meat_hac(scores, weights, m, status, adjust = .false.)
   call assert_true(status == SANDWICH_SUCCESS, 'HAC status')
   call assert_close(m(1, 1), 0.8_dp, 1.0e-13_dp, 'HAC 11')
   call assert_close(m(1, 2), -0.1_dp, 1.0e-13_dp, 'HAC 12')
   call assert_close(m(2, 2), 2.0_dp, 1.0e-13_dp, 'HAC 22')

   do i = 1, 40
      series(i, 1) = sin(0.21_dp * real(i, dp)) + 0.02_dp * real(i, dp)
      series(i, 2) = cos(0.13_dp * real(i, dp)) - 0.01_dp * real(i, dp)
   end do
   call prewhite_var(series, 1, residuals, recolor, status)
   call assert_true(status == SANDWICH_SUCCESS, 'prewhite status')
   call assert_true(all(shape(residuals) == [39, 2]), 'prewhite residual shape')
   call assert_true(all(shape(recolor) == [2, 2]), 'prewhite recolor shape')

   call bandwidth_andrews(series, 'Quadratic Spectral', bw, status, prewhite_order = 0)
   call assert_true(status == SANDWICH_SUCCESS .and. bw > 0.0_dp, 'Andrews bandwidth')
   call bandwidth_newey_west(series, 'Bartlett', bw, status, prewhite_order = 0)
   call assert_true(status == SANDWICH_SUCCESS .and. bw >= 0.0_dp, 'Newey-West bandwidth')

   print '(a)', 'test_kernels_hac: PASS'

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         print '(a)', 'FAIL: ' // trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: message
      if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
         print '(a,2(1x,es24.16))', 'FAIL: ' // trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_vector_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      character(len=*), intent(in) :: message
      integer :: j
      call assert_true(size(actual) == size(expected), trim(message) // ' size')
      do j = 1, size(actual)
         call assert_close(actual(j), expected(j), tolerance, trim(message))
      end do
   end subroutine assert_vector_close

end program test_kernels_hac
