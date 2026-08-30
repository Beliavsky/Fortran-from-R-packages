! SPDX-License-Identifier: GPL-2.0-or-later
program advanced_api_example
   use gbm3
   implicit none

   integer, parameter :: n = 60
   real(dp) :: x(n, 2), y(n)
   real(dp), allocatable :: staged(:, :), importance(:)
   integer :: counts(3), i
   type(gbm_options) :: options
   type(gbm_model) :: model
   type(gbm_cv_result) :: cv

   do i = 1, n
      x(i, 1) = -2.0_dp + 4.0_dp * real(i - 1, dp) / real(n - 1, dp)
      x(i, 2) = sin(0.25_dp * real(i, dp))
      y(i) = 1.5_dp * x(i, 1) - 0.4_dp * x(i, 2)
   end do

   options = gbm_options(distribution=GBM_GAUSSIAN, num_trees=20, &
                         interaction_depth=2, min_num_obs_in_node=3, &
                         shrinkage=0.08_dp, bag_fraction=0.7_dp)
   call gbm_set_seed(1234)
   call gbm_fit(x, y, model, options)
   call gbm_continue(model, x, y, 10)

   counts = [10, 20, 30]
   staged = gbm_predict_staged(model, x, counts)
   call gbm_set_seed(5678)
   call gbm_permutation_importance(model, x, y, importance, rescale=.true.)

   options%num_trees = 20
   call gbm_set_seed(2468)
   call gbm_cross_validate(x, y, cv, 3, options)

   print '(a,i0)', 'Trees after continuation: ', model%n_trees
   print '(a,3f12.6)', 'First staged prediction: ', staged(1, :)
   print '(a,2f12.6)', 'Permutation importance: ', importance
   print '(a,i0)', 'Best CV iteration: ', cv%best_iteration
end program advanced_api_example
