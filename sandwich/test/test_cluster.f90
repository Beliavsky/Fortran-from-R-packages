! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_cluster
   use sandwich, only : dp, SANDWICH_SUCCESS, meat_cluster, ols_model, fit_ols
   implicit none

   real(dp) :: scores(6, 2), x(12, 2), y(12)
   integer :: cluster1(6, 1), cluster2(6, 2), cluster_hc(12, 1), status, i
   real(dp), allocatable :: m(:, :)
   type(ols_model) :: model

   scores = reshape([ &
      1.0_dp, -1.0_dp, 2.0_dp, 0.0_dp, 1.0_dp, -2.0_dp, &
      2.0_dp,  1.0_dp, 0.0_dp,-2.0_dp, 1.0_dp,  1.0_dp], [6, 2])
   cluster1(:, 1) = [1, 1, 2, 2, 3, 3]
   cluster2(:, 1) = cluster1(:, 1)
   cluster2(:, 2) = [1, 2, 1, 2, 1, 2]

   call meat_cluster(scores, cluster1, m, status, type = 'HC0', cadjust = .false.)
   call assert_true(status == SANDWICH_SUCCESS, 'one-way cluster status')
   call assert_close(m(1, 1), 5.0_dp / 6.0_dp, 1.0e-13_dp, 'one-way 11')
   call assert_close(m(1, 2), -1.0_dp, 1.0e-13_dp, 'one-way 12')
   call assert_close(m(2, 2), 17.0_dp / 6.0_dp, 1.0e-13_dp, 'one-way 22')

   call meat_cluster(scores, cluster2, m, status, type = 'HC0', cadjust = .false.)
   call assert_true(status == SANDWICH_SUCCESS, 'two-way cluster status')
   call assert_close(m(1, 1), 19.0_dp / 6.0_dp, 1.0e-13_dp, 'two-way 11')
   call assert_close(m(1, 2), 1.0_dp, 1.0e-13_dp, 'two-way 12')
   call assert_close(m(2, 2), 2.5_dp, 1.0e-13_dp, 'two-way 22')

   call meat_cluster(scores, cluster2, m, status, type = 'HC0', cadjust = .false., multi0 = .true.)
   call assert_true(status == SANDWICH_SUCCESS, 'multi0 status')
   call assert_close(m(1, 1), 19.0_dp / 6.0_dp, 1.0e-13_dp, 'multi0 11')

   do i = 1, 12
      x(i, :) = [1.0_dp, real(i - 1, dp)]
      y(i) = 1.0_dp + 0.4_dp * x(i, 2) + 0.2_dp * sin(real(i, dp))
      cluster_hc(i, 1) = (i + 2) / 3
   end do
   call fit_ols(x, y, model, status)
   call assert_true(status == SANDWICH_SUCCESS, 'cluster HC fit')
   call meat_cluster(model%scores, cluster_hc, m, status, type = 'HC2', &
      x = x, residuals = model%residuals, hat = model%hat)
   call assert_true(status == SANDWICH_SUCCESS, 'cluster HC2 status')
   call assert_true(maxval(abs(m - transpose(m))) < 1.0e-12_dp, 'cluster HC2 symmetry')
   call meat_cluster(model%scores, cluster_hc, m, status, type = 'HC3', &
      x = x, residuals = model%residuals, hat = model%hat)
   call assert_true(status == SANDWICH_SUCCESS, 'cluster HC3 status')
   call assert_true(maxval(abs(m - transpose(m))) < 1.0e-12_dp, 'cluster HC3 symmetry')

   print '(a)', 'test_cluster: PASS'

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

end program test_cluster
