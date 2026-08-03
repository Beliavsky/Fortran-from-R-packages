! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_panel
   use sandwich, only : dp, SANDWICH_SUCCESS, meat_panel_longitudinal, meat_panel_corrected
   implicit none

   real(dp) :: scores(12, 2), x(12, 2), residuals(12)
   real(dp), parameter :: time_effect(4) = [0.2_dp, -0.1_dp, 0.3_dp, 0.0_dp]
   real(dp), parameter :: cluster_effect(3) = [0.05_dp, -0.03_dp, 0.02_dp]
   integer :: cluster(12), time(12), status, i, c, t
   real(dp), allocatable :: m(:, :)

   i = 0
   do c = 1, 3
      do t = 1, 4
         i = i + 1
         cluster(i) = c
         time(i) = t
         scores(i, 1) = 0.2_dp * real(c, dp) + 0.1_dp * real(t, dp)
         scores(i, 2) = (-1.0_dp)**(c + t) * &
            (0.3_dp * real(c, dp) - 0.05_dp * real(t, dp))
      end do
   end do

   call meat_panel_longitudinal(scores, cluster, time, m, status, kernel = 'Bartlett', &
      lag = 1, bandwidth = 2.0_dp, adjust = .false., aggregate = .true.)
   call assert_true(status == SANDWICH_SUCCESS, 'PL aggregate status')
   call assert_close(m(1, 1), 2.265_dp, 1.0e-13_dp, 'PL aggregate 11')
   call assert_close(m(1, 2), 0.0025_dp, 1.0e-13_dp, 'PL aggregate 12')
   call assert_close(m(2, 2), 0.0195833333333333_dp, 1.0e-13_dp, 'PL aggregate 22')

   call meat_panel_longitudinal(scores, cluster, time, m, status, kernel = 'Bartlett', &
      lag = 1, bandwidth = 2.0_dp, adjust = .false., aggregate = .false.)
   call assert_true(status == SANDWICH_SUCCESS, 'PL nonaggregate status')
   call assert_close(m(1, 1), 0.801666666666667_dp, 1.0e-13_dp, 'PL nonaggregate 11')
   call assert_close(m(1, 2), 0.000833333333333333_dp, 1.0e-13_dp, 'PL nonaggregate 12')
   call assert_close(m(2, 2), 0.07375_dp, 1.0e-13_dp, 'PL nonaggregate 22')

   i = 0
   do t = 1, 4
      do c = 1, 3
         i = i + 1
         cluster(i) = c
         time(i) = t
         x(i, 1) = 1.0_dp
         x(i, 2) = 0.2_dp * real(c, dp) + 0.1_dp * real(t, dp)
         residuals(i) = time_effect(t) + cluster_effect(c)
      end do
   end do
   call meat_panel_corrected(x, residuals, cluster, time, m, status)
   call assert_true(status == SANDWICH_SUCCESS, 'PC status')
   call assert_close(m(1, 1), 0.113533333333333_dp, 1.0e-13_dp, 'PC 11')
   call assert_close(m(1, 2), 0.0731166666666667_dp, 1.0e-13_dp, 'PC 12')
   call assert_close(m(2, 2), 0.048515_dp, 1.0e-13_dp, 'PC 22')

   print '(a)', 'test_panel: PASS'

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

end program test_panel
