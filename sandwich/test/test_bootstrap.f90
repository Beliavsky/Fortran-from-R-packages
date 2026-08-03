! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_bootstrap
   use sandwich, only : dp, SANDWICH_SUCCESS, bootstrap_covariance, jackknife_covariance, &
      vcov_bootstrap_ols
   implicit none

   real(dp) :: replicates(4, 2), x(8, 2), y(8)
   integer :: cluster(8, 1), status, i
   real(dp), allocatable :: cov(:, :), cov2(:, :)

   replicates = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
                         2.0_dp, 4.0_dp, 3.0_dp, 5.0_dp], [4, 2])
   call bootstrap_covariance(replicates, cov, status)
   call assert_true(status == SANDWICH_SUCCESS, 'bootstrap covariance status')
   call assert_close(cov(1, 1), 5.0_dp / 3.0_dp, 1.0e-13_dp, 'bootstrap 11')
   call assert_close(cov(1, 2), 4.0_dp / 3.0_dp, 1.0e-13_dp, 'bootstrap 12')
   call assert_close(cov(2, 2), 5.0_dp / 3.0_dp, 1.0e-13_dp, 'bootstrap 22')

   call jackknife_covariance(replicates, cov, status)
   call assert_true(status == SANDWICH_SUCCESS, 'jackknife covariance status')
   call assert_close(cov(1, 1), 3.75_dp, 1.0e-13_dp, 'jackknife 11')
   call assert_close(cov(1, 2), 3.0_dp, 1.0e-13_dp, 'jackknife 12')
   call assert_close(cov(2, 2), 3.75_dp, 1.0e-13_dp, 'jackknife 22')

   do i = 1, 8
      x(i, 1) = 1.0_dp
      x(i, 2) = real(i - 1, dp)
      y(i) = 1.0_dp + 0.7_dp * x(i, 2) + 0.2_dp * sin(real(i, dp))
      cluster(i, 1) = (i + 1) / 2
   end do
   call vcov_bootstrap_ols(x, y, cluster, cov, status, replications = 40, &
      type = 'rademacher', seed = 12345, fix = .true.)
   call assert_true(status == SANDWICH_SUCCESS, 'wild bootstrap status')
   call assert_true(all(shape(cov) == [2, 2]), 'wild bootstrap shape')
   call assert_true(maxval(abs(cov - transpose(cov))) < 1.0e-12_dp, 'wild bootstrap symmetry')
   call assert_true(all([(cov(i, i) >= 0.0_dp, i = 1, 2)]), 'wild bootstrap diagonal')

   call vcov_bootstrap_ols(x, y, cluster, cov2, status, replications = 40, &
      type = 'rademacher', seed = 12345, fix = .true.)
   call assert_true(status == SANDWICH_SUCCESS, 'repeat wild bootstrap status')
   call assert_true(maxval(abs(cov - cov2)) < 1.0e-14_dp, 'bootstrap reproducibility')

   print '(a)', 'test_bootstrap: PASS'

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

end program test_bootstrap
