! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program cluster_panel_example
   use sandwich, only : dp, ols_model, fit_ols, vcov_cluster, &
      vcov_panel_longitudinal, SANDWICH_SUCCESS
   implicit none

   integer, parameter :: firms = 5, periods = 6, n = firms * periods
   real(dp) :: x(n, 2), y(n)
   integer :: cluster_matrix(n, 1), cluster(n), time(n), i, c, t, status
   real(dp), allocatable :: clustered(:, :), panel(:, :)
   type(ols_model) :: model

   i = 0
   do c = 1, firms
      do t = 1, periods
         i = i + 1
         cluster(i) = c
         cluster_matrix(i, 1) = c
         time(i) = t
         x(i, :) = [1.0_dp, 0.15_dp * real(c, dp) + 0.1_dp * real(t, dp)]
         y(i) = 0.7_dp + 1.2_dp * x(i, 2) + 0.1_dp * sin(real(c + t, dp))
      end do
   end do
   call fit_ols(x, y, model, status)
   if (status /= SANDWICH_SUCCESS) error stop
   call vcov_cluster(model%scores, cluster_matrix, model%bread, clustered, status)
   if (status /= SANDWICH_SUCCESS) error stop
   call vcov_panel_longitudinal(model%scores, cluster, time, model%bread, panel, status, &
      lag = 2, adjust = .false.)
   if (status /= SANDWICH_SUCCESS) error stop

   print '(a,2f12.6)', 'cluster variances: ', clustered(1, 1), clustered(2, 2)
   print '(a,2f12.6)', 'panel variances:   ', panel(1, 1), panel(2, 2)
end program cluster_panel_example
