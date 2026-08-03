! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program bootstrap_example
   use sandwich, only : dp, vcov_bootstrap_ols, SANDWICH_SUCCESS
   implicit none

   real(dp) :: x(20, 2), y(20)
   integer :: cluster(20, 1), i, status
   real(dp), allocatable :: covariance(:, :)

   do i = 1, 20
      x(i, :) = [1.0_dp, real(i - 1, dp) / 5.0_dp]
      y(i) = 1.0_dp + 0.4_dp * x(i, 2) + 0.15_dp * sin(real(i, dp))
      cluster(i, 1) = (i + 3) / 4
   end do

   call vcov_bootstrap_ols(x, y, cluster, covariance, status, replications = 200, &
      type = 'rademacher', seed = 20260801, fix = .true.)
   if (status /= SANDWICH_SUCCESS) error stop
   print '(a,2f12.6)', 'wild-bootstrap variances: ', covariance(1, 1), covariance(2, 2)
end program bootstrap_example
