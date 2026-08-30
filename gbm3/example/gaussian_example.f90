! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of gbm3 3.0.3; see NOTICE.md and upstream/.
program gaussian_example
   use gbm3
   implicit none
   integer, parameter :: n = 100
   real(dp) :: x(n, 2), y(n), pred(n)
   real(dp), allocatable :: influence(:)
   type(gbm_options) :: options
   type(gbm_model) :: model
   integer :: i

   do i = 1, n
      x(i, 1) = -2.0_dp + 4.0_dp * real(i - 1, dp) / real(n - 1, dp)
      x(i, 2) = cos(0.2_dp * real(i, dp))
      y(i) = 1.5_dp + 2.0_dp * x(i, 1) - 0.3_dp * x(i, 2)
   end do

   options = gbm_options(distribution=GBM_GAUSSIAN, num_trees=100, interaction_depth=2, &
                         min_num_obs_in_node=5, shrinkage=0.05_dp, bag_fraction=0.5_dp)
   call gbm_set_seed(1234)
   call gbm_fit(x, y, model, options)
   pred = gbm_predict(model, x)
   call gbm_relative_influence(model, influence, normalize=.true.)

   print '(a,f12.6)', 'Initial estimate: ', model%init_f
   print '(a,f12.6)', 'Final training error: ', model%train_error(model%n_trees)
   print '(a,f12.6)', 'Prediction MSE: ', sum((pred - y) ** 2) / real(n, dp)
   print '(a,2f10.3)', 'Relative influence (%): ', influence
end program gaussian_example
