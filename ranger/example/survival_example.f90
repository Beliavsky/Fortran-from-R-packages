! SPDX-License-Identifier: GPL-3.0-only
! See NOTICE.md for ranger upstream provenance.
program survival_example
   use r_kinds, only : dp
   use ranger, only : ranger_options, ranger_survival_forest
   use ranger, only : fit_ranger_survival, predict_ranger_survival
   implicit none
   real(dp) :: x(20, 1), time(20), grid(6), chf(20, 6), survival(20, 6)
   integer :: event(20), i, status
   type(ranger_options) :: options
   type(ranger_survival_forest) :: forest

   do i = 1, 10
      x(i, 1) = 0.0_dp
      time(i) = 1.0_dp + 0.25_dp * real(i - 1, dp)
      event(i) = 1
   end do
   do i = 11, 20
      x(i, 1) = 1.0_dp
      time(i) = 5.0_dp + 0.25_dp * real(i - 11, dp)
      event(i) = merge(1, 0, mod(i, 3) /= 0)
   end do
   grid = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]

   options%num_trees = 25
   options%mtry = 1
   options%min_node_size = 3
   options%seed = 2026
   call fit_ranger_survival(x, time, event, forest, options=options, time_interest=grid, status=status)
   if (status /= 0) error stop 'survival fit failed'
   call predict_ranger_survival(forest, x, chf, survival)

   print '(a,6f9.4)', 'survival for first case: ', survival(1, :)
end program survival_example
