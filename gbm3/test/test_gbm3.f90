! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of gbm3 3.0.3; see NOTICE.md and upstream/.
program test_gbm3
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
   use gbm3
   implicit none

   call test_distribution_initial_values()
   call test_gaussian_fit()
   call test_bernoulli_fit()
   call test_other_simple_distributions()
   call test_categorical_and_missing_tree()
   call test_pairwise_fit()
   call test_cox_fit()
   call test_counting_cox_and_baseline()
   call test_validation_split()
   call test_continue_and_staged_prediction()
   call test_cross_validation()
   call test_diagnostics_api()
   call test_pairwise_and_cox_expanded_apis()
   print '(a)', 'gbm3 tests passed'

contains

   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(*), intent(in) :: message
      if (.not. condition) then
         print '(a)', 'FAIL: ' // trim(message)
         error stop 1
      end if
   end subroutine check

   subroutine test_distribution_initial_values()
      use gbm3_distributions, only : dist_init_f
      real(dp) :: y(4), off(4), w(4), v
      type(gbm_options) :: opt
      y = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      off = 0.0_dp
      w = 1.0_dp

      opt = gbm_options(distribution=GBM_GAUSSIAN)
      v = dist_init_f(y, off, w, opt)
      call check(abs(v - 2.5_dp) < 1.0e-12_dp, 'Gaussian initial estimate')

      opt%distribution = GBM_POISSON
      v = dist_init_f(y, off, w, opt)
      call check(abs(v - log(2.5_dp)) < 1.0e-12_dp, 'Poisson initial estimate')

      opt%distribution = GBM_GAMMA
      v = dist_init_f(y, off, w, opt)
      call check(abs(v - log(2.5_dp)) < 1.0e-12_dp, 'Gamma initial estimate')

      y = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
      opt%distribution = GBM_BERNOULLI
      v = dist_init_f(y, off, w, opt)
      call check(abs(v) < 1.0e-12_dp, 'Bernoulli initial estimate')

      opt%distribution = GBM_ADABOOST
      v = dist_init_f(y, off, w, opt)
      call check(abs(v) < 1.0e-12_dp, 'AdaBoost initial estimate')
   end subroutine test_distribution_initial_values

   subroutine test_gaussian_fit()
      integer, parameter :: n = 80
      real(dp) :: x(n, 2), y(n), pred(n), grid(5, 1)
      real(dp), allocatable :: influence(:), pd(:)
      type(gbm_options) :: opt
      type(gbm_model) :: model
      integer :: i

      do i = 1, n
         x(i, 1) = -2.0_dp + 4.0_dp * real(i - 1, dp) / real(n - 1, dp)
         x(i, 2) = sin(0.3_dp * real(i, dp))
         y(i) = 2.0_dp * x(i, 1) - 0.5_dp * x(i, 2) + 0.15_dp * cos(real(i, dp))
      end do
      opt = gbm_options(distribution=GBM_GAUSSIAN, num_trees=40, interaction_depth=2, &
                        min_num_obs_in_node=3, shrinkage=0.1_dp, bag_fraction=0.65_dp)
      call gbm_set_seed(12345)
      call gbm_fit(x, y, model, opt)
      pred = gbm_predict(model, x)
      call check(model%n_trees == 40, 'Gaussian tree count')
      call check(model%train_error(40) < model%train_error(1), 'Gaussian deviance improves')
      call check(sum((pred - y) ** 2) / real(n, dp) < 0.35_dp, 'Gaussian prediction MSE')
      call gbm_relative_influence(model, influence, normalize=.true.)
      call check(abs(sum(influence) - 100.0_dp) < 1.0e-9_dp, 'relative influence normalization')
      call check(influence(1) > influence(2), 'dominant Gaussian feature has greater influence')
      grid(:, 1) = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
      pd = gbm_partial_dependence(model, grid, [1])
      call check(size(pd) == 5 .and. all(ieee_is_finite(pd)), 'partial dependence')
      call check(pd(5) > pd(1), 'partial dependence direction')
      call check(all(ieee_is_nan(model%validation_error)), 'no-validation error is NaN')
   end subroutine test_gaussian_fit

   subroutine test_bernoulli_fit()
      integer, parameter :: n = 100
      real(dp) :: x(n, 1), y(n), p(n)
      type(gbm_options) :: opt
      type(gbm_model) :: model
      integer :: i

      do i = 1, n
         x(i, 1) = -3.0_dp + 6.0_dp * real(i - 1, dp) / real(n - 1, dp)
         if (x(i, 1) > 0.15_dp) then
            y(i) = 1.0_dp
         else
            y(i) = 0.0_dp
         end if
      end do
      opt = gbm_options(distribution=GBM_BERNOULLI, num_trees=30, interaction_depth=1, &
                        min_num_obs_in_node=3, shrinkage=0.15_dp, bag_fraction=0.7_dp)
      call gbm_set_seed(222)
      call gbm_fit(x, y, model, opt)
      p = gbm_predict_response(model, x)
      call check(model%train_error(30) < model%train_error(1), 'Bernoulli deviance improves')
      call check(minval(p) >= 0.0_dp .and. maxval(p) <= 1.0_dp, 'Bernoulli response probabilities')
      call check(sum(merge(1.0_dp, 0.0_dp, (p >= 0.5_dp) .eqv. (y > 0.5_dp))) / real(n, dp) > 0.9_dp, &
                 'Bernoulli classification accuracy')
   end subroutine test_bernoulli_fit

   subroutine test_other_simple_distributions()
      integer, parameter :: n = 48
      real(dp) :: x(n, 1), y(n), pred(n)
      type(gbm_options) :: opt
      type(gbm_model) :: model
      integer :: i, dist
      integer, parameter :: dists(8) = [GBM_POISSON, GBM_GAMMA, GBM_LAPLACE, GBM_TDIST, &
                                        GBM_QUANTILE, GBM_ADABOOST, GBM_HUBERIZED, GBM_TWEEDIE]

      do dist = 1, size(dists)
         do i = 1, n
            x(i, 1) = -1.5_dp + 3.0_dp * real(i - 1, dp) / real(n - 1, dp)
         end do
         select case (dists(dist))
         case (GBM_POISSON)
            y = max(0.0_dp, anint(exp(0.4_dp + 0.45_dp * x(:, 1))))
         case (GBM_GAMMA, GBM_TWEEDIE)
            y = exp(0.2_dp + 0.35_dp * x(:, 1)) + 0.1_dp
         case (GBM_ADABOOST, GBM_HUBERIZED)
            y = merge(1.0_dp, 0.0_dp, x(:, 1) > 0.0_dp)
         case default
            y = 1.4_dp * x(:, 1) + 0.08_dp * sin([(real(i, dp), i=1,n)])
         end select
         opt = gbm_options(distribution=dists(dist), num_trees=8, interaction_depth=1, &
                           min_num_obs_in_node=2, shrinkage=0.08_dp, bag_fraction=0.65_dp)
         call gbm_set_seed(1000 + dist)
         call gbm_fit(x, y, model, opt)
         pred = gbm_predict(model, x)
         call check(all(ieee_is_finite(pred)), 'finite simple-distribution predictions')
         call check(ieee_is_finite(model%train_error(model%n_trees)), 'finite simple-distribution deviance')
      end do
   end subroutine test_other_simple_distributions

   subroutine test_categorical_and_missing_tree()
      use gbm3_tree, only : grow_tree, adjust_tree
      use gbm3_types, only : gbm_tree
      integer, parameter :: n = 18
      real(dp) :: x(n, 1), r(n), w(n), delta(n)
      logical :: bag(n)
      integer :: vc(1), mono(1), i
      integer, allocatable :: assignment(:)
      type(gbm_tree) :: tree

      do i = 1, n
         x(i, 1) = real(mod(i - 1, 3), dp)
         r(i) = merge(2.0_dp, -1.0_dp, mod(i - 1, 3) == 2)
      end do
      x(3, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
      x(12, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
      w = 1.0_dp
      bag = .true.
      vc = 3
      mono = 0
      call gbm_set_seed(42)
      call grow_tree(x, r, w, bag, vc, mono, 1, 2, 1, 0.1_dp, tree, assignment)
      call check(.not. tree%nodes(1)%is_terminal, 'categorical predictor can split')
      call check(allocated(tree%nodes(1)%left_categories), 'categorical split stores categories')
      call adjust_tree(tree, 2, assignment, delta)
      call check(all(ieee_is_finite(delta)), 'categorical/missing adjustment is finite')
   end subroutine test_categorical_and_missing_tree

   subroutine test_pairwise_fit()
      integer, parameter :: ng = 10, ni = 4, n = ng * ni
      real(dp) :: x(n, 2), y(n), pred(n)
      integer :: group(n), g, k, i
      type(gbm_options) :: opt
      type(gbm_model) :: model

      i = 0
      do g = 1, ng
         do k = 1, ni
            i = i + 1
            group(i) = g
            y(i) = real(ni - k, dp)
            x(i, 1) = y(i) + 0.02_dp * real(g, dp)
            x(i, 2) = real(mod(g + k, 3), dp)
         end do
      end do
      opt = gbm_options(distribution=GBM_PAIRWISE, num_trees=15, interaction_depth=1, &
                        min_num_obs_in_node=2, shrinkage=0.08_dp, bag_fraction=0.6_dp, &
                        pairwise_metric=GBM_METRIC_NDCG)
      call gbm_set_seed(77)
      call gbm_fit(x, y, model, opt, group=group)
      pred = gbm_predict(model, x)
      call check(all(ieee_is_finite(pred)), 'pairwise predictions are finite')
      call check(ieee_is_finite(model%train_error(model%n_trees)), 'pairwise deviance is finite')
      call check(model%train_error(model%n_trees) <= model%train_error(1) + 1.0e-10_dp, &
                 'pairwise ranking error does not worsen')
   end subroutine test_pairwise_fit

   subroutine test_cox_fit()
      integer, parameter :: n = 36
      real(dp) :: x(n, 1), surv(n, 2), pred(n)
      integer :: strata(n), i
      type(gbm_options) :: opt
      type(gbm_model) :: model

      do i = 1, n
         x(i, 1) = -1.5_dp + 3.0_dp * real(i - 1, dp) / real(n - 1, dp)
         surv(i, 1) = 8.0_dp - 1.5_dp * x(i, 1) + 0.04_dp * real(mod(i, 5), dp)
         surv(i, 2) = merge(1.0_dp, 0.0_dp, mod(i, 4) /= 0)
         strata(i) = 1
      end do
      opt = gbm_options(distribution=GBM_COXPH, num_trees=8, interaction_depth=1, &
                        min_num_obs_in_node=3, shrinkage=0.05_dp, bag_fraction=0.7_dp)
      call gbm_set_seed(404)
      call gbm_fit(x, surv, model, opt, strata=strata)
      pred = gbm_predict(model, x)
      call check(all(ieee_is_finite(pred)), 'Cox predictions are finite')
      call check(all(ieee_is_finite(model%train_error)), 'Cox training deviance is finite')
   end subroutine test_cox_fit

   subroutine test_counting_cox_and_baseline()
      integer, parameter :: n = 24
      real(dp) :: x(n, 1), surv(n, 3), pred(n), weight(n), eta(n)
      real(dp), allocatable :: times(:), hazard(:), cumulative(:)
      integer :: strata(n), i
      type(gbm_options) :: opt
      type(gbm_model) :: model

      do i = 1, n
         x(i, 1) = real(i - 1, dp) / real(n - 1, dp)
         surv(i, 1) = 0.05_dp * real(mod(i - 1, 3), dp)
         surv(i, 2) = 1.0_dp + 0.1_dp * real(i, dp)
         surv(i, 3) = merge(1.0_dp, 0.0_dp, mod(i, 4) /= 0)
         strata(i) = 1
      end do
      weight = 1.0_dp
      opt = gbm_options(distribution=GBM_COXPH, num_trees=4, interaction_depth=1, &
                        min_num_obs_in_node=2, shrinkage=0.05_dp, bag_fraction=0.7_dp)
      call gbm_set_seed(505)
      call gbm_fit(x, surv, model, opt, strata=strata)
      pred = gbm_predict(model, x)
      call check(all(ieee_is_finite(pred)), 'counting-process Cox predictions')
      eta = pred
      call cox_baseline_hazard(surv, strata, eta, weight, GBM_TIES_EFRON, times, hazard, cumulative)
      call check(size(times) > 0, 'Cox baseline hazard has event times')
      call check(all(hazard >= 0.0_dp) .and. all(ieee_is_finite(hazard)), 'Cox baseline hazards are finite')
      call check(all(ieee_is_finite(cumulative)), 'Cox cumulative baseline hazard is finite')
   end subroutine test_counting_cox_and_baseline

   subroutine test_validation_split()
      integer, parameter :: n = 50
      real(dp) :: x(n, 1), y(n)
      type(gbm_options) :: opt
      type(gbm_model) :: model
      integer :: i, best

      do i = 1, n
         x(i, 1) = real(i, dp) / real(n, dp)
         y(i) = 3.0_dp * x(i, 1) - 0.5_dp
      end do
      opt = gbm_options(distribution=GBM_GAUSSIAN, num_trees=6, interaction_depth=1, &
                        min_num_obs_in_node=2, shrinkage=0.1_dp, bag_fraction=0.7_dp, num_train=35)
      call gbm_set_seed(606)
      call gbm_fit(x, y, model, opt)
      call check(all(ieee_is_finite(model%validation_error)), 'validation errors are computed')
      best = gbm_best_iteration(model, 'test')
      call check(best >= 1 .and. best <= model%n_trees, 'validation best-iteration selector')
   end subroutine test_validation_split


   subroutine test_continue_and_staged_prediction()
      integer, parameter :: n = 72
      real(dp) :: x(n, 2), y(n)
      real(dp), allocatable :: p_full(:), p_more(:), p_stage(:, :), contrib(:, :), p_one(:)
      type(gbm_options) :: opt_full, opt_part
      type(gbm_model) :: full_model, more_model
      integer, parameter :: counts(4) = [0, 1, 10, 40]
      integer :: i, j

      do i = 1, n
         x(i, 1) = -1.8_dp + 3.6_dp * real(i - 1, dp) / real(n - 1, dp)
         x(i, 2) = cos(0.17_dp * real(i, dp))
         y(i) = 1.3_dp * x(i, 1) - 0.7_dp * x(i, 2) + 0.05_dp * sin(real(i, dp))
      end do
      opt_full = gbm_options(distribution=GBM_GAUSSIAN, num_trees=40, interaction_depth=2, &
                             min_num_obs_in_node=3, shrinkage=0.08_dp, bag_fraction=0.65_dp)
      opt_part = opt_full
      opt_part%num_trees = 20

      call gbm_set_seed(9001)
      call gbm_fit(x, y, full_model, opt_full)
      call gbm_set_seed(9001)
      call gbm_fit(x, y, more_model, opt_part)
      call gbm_continue(more_model, x, y, 20)

      p_full = gbm_predict(full_model, x)
      p_more = gbm_predict(more_model, x)
      call check(more_model%n_trees == 40 .and. more_model%options%num_trees == 40, &
                 'continuation updates tree counts')
      call check(more_model%n_rows == n .and. more_model%n_train == n, 'continuation retains fit dimensions')
      call check(maxval(abs(p_full - p_more)) < 1.0e-12_dp, 'continuation matches one-shot fit')
      call check(maxval(abs(full_model%train_error - more_model%train_error)) < 1.0e-12_dp, &
                 'continuation matches one-shot training history')

      p_stage = gbm_predict_staged(more_model, x, counts)
      do j = 1, size(counts)
         p_one = gbm_predict(more_model, x, n_trees=counts(j))
         call check(maxval(abs(p_stage(:, j) - p_one)) < 1.0e-13_dp, 'staged prediction matches scalar tree count')
      end do
      contrib = gbm_predict_trees(more_model, x)
      call check(size(contrib, 2) == more_model%n_trees, 'tree contribution matrix tree count')
      call check(maxval(abs(more_model%init_f + sum(contrib, dim=2) - p_more)) < 1.0e-13_dp, &
                 'tree contributions reconstruct link prediction')
   end subroutine test_continue_and_staged_prediction


   subroutine test_cross_validation()
      integer, parameter :: n = 48
      real(dp) :: x(n, 2), y(n), yb(n)
      integer :: id(n), i, f
      type(gbm_options) :: opt
      type(gbm_cv_result) :: cv, cvb
      type(gbm_model) :: full_model

      do i = 1, n
         x(i, 1) = -2.0_dp + 4.0_dp * real(i - 1, dp) / real(n - 1, dp)
         x(i, 2) = sin(0.21_dp * real(i, dp))
         y(i) = 1.7_dp * x(i, 1) + 0.3_dp * x(i, 2)
         id(i) = (i + 1) / 2
         yb(i) = merge(1.0_dp, 0.0_dp, mod(id(i), 2) == 0)
      end do
      opt = gbm_options(distribution=GBM_GAUSSIAN, num_trees=12, interaction_depth=1, &
                        min_num_obs_in_node=2, shrinkage=0.1_dp, bag_fraction=0.7_dp)
      call gbm_set_seed(3333)
      call gbm_cross_validate(x, y, cv, 4, opt, id=id, full_model=full_model)
      call check(cv%n_folds == 4 .and. size(cv%error) == 12, 'cross-validation dimensions')
      call check(cv%best_iteration >= 1 .and. cv%best_iteration <= 12, 'cross-validation best iteration')
      call check(gbm_cv_best_iteration(cv) == cv%best_iteration, 'cross-validation best iteration accessor')
      call check(all(ieee_is_finite(cv%error)) .and. all(ieee_is_finite(cv%fitted)), &
                 'cross-validation finite errors and fitted values')
      call check(full_model%n_trees == 12, 'cross-validation returns full-data model')
      do i = 2, n, 2
         call check(cv%fold_id(i - 1) == cv%fold_id(i), 'repeated ids remain in one CV fold')
      end do

      opt%distribution = GBM_BERNOULLI
      opt%num_trees = 8
      call gbm_set_seed(4444)
      call gbm_cross_validate(x, yb, cvb, 4, opt, id=id, stratify=.true.)
      do f = 1, 4
         call check(any((cvb%fold_id == f) .and. (yb > 0.5_dp)), 'stratified fold contains positives')
         call check(any((cvb%fold_id == f) .and. (yb < 0.5_dp)), 'stratified fold contains negatives')
      end do
   end subroutine test_cross_validation


   subroutine test_diagnostics_api()
      integer, parameter :: n = 80
      real(dp) :: x(n, 2), y(n)
      real(dp), allocatable :: importance(:)
      type(gbm_options) :: opt
      type(gbm_model) :: model
      type(gbm_tree) :: tree
      type(gbm_node) :: node
      integer :: i, best
      real(dp) :: h

      do i = 1, n
         x(i, 1) = -2.0_dp + 4.0_dp * real(i - 1, dp) / real(n - 1, dp)
         x(i, 2) = sin(0.61_dp * real(i, dp))
         y(i) = 2.2_dp * x(i, 1) + 0.02_dp * cos(real(i, dp))
      end do
      opt = gbm_options(distribution=GBM_GAUSSIAN, num_trees=25, interaction_depth=2, &
                        min_num_obs_in_node=3, shrinkage=0.1_dp, bag_fraction=0.7_dp)
      call gbm_set_seed(8181)
      call gbm_fit(x, y, model, opt)
      best = gbm_best_iteration(model, 'train')
      call check(best >= 1 .and. best <= model%n_trees, 'train best-iteration selector')
      best = gbm_best_iteration(model, 'oob_raw')
      call check(best >= 1 .and. best <= model%n_trees, 'raw OOB best-iteration selector')

      call gbm_set_seed(9191)
      call gbm_permutation_importance(model, x, y, importance)
      call check(size(importance) == 2 .and. all(ieee_is_finite(importance)), 'permutation importance dimensions')
      call check(importance(1) > importance(2), 'permutation importance identifies dominant feature')
      call gbm_set_seed(9191)
      call gbm_permutation_importance(model, x, y, importance, rescale=.true.)
      call check(abs(maxval(importance) - 1.0_dp) < 1.0e-12_dp, 'permutation importance rescaling')

      tree = gbm_get_tree(model, 1)
      call check(tree%n_nodes == gbm_num_nodes(model, 1), 'tree inspection node count')
      node = gbm_get_node(model, 1, 1)
      call check(node%num_obs == tree%nodes(1)%num_obs, 'node inspection accessor')
      h = gbm_interaction_strength(model, x, [1, 2])
      call check(ieee_is_finite(h) .and. h >= 0.0_dp .and. h <= 1.0_dp, 'Friedman H interaction strength')
   end subroutine test_diagnostics_api


   subroutine test_pairwise_and_cox_expanded_apis()
      integer, parameter :: ng = 8, ni = 3, nr = ng * ni, ns = 27
      real(dp) :: xr(nr, 2), yr(nr), pr_full(nr), pr_more(nr)
      real(dp) :: xs(ns, 1), surv(ns, 2), ps_full(ns), ps_more(ns)
      integer :: group(nr), strata(ns), g, k, i
      type(gbm_options) :: opt, part
      type(gbm_model) :: pair_full, pair_more, cox_full, cox_more
      type(gbm_cv_result) :: cv_pair, cv_cox

      i = 0
      do g = 1, ng
         do k = 1, ni
            i = i + 1
            group(i) = g
            yr(i) = real(ni - k, dp)
            xr(i, 1) = yr(i) + 0.01_dp * real(g, dp)
            xr(i, 2) = real(mod(g + 2 * k, 4), dp)
         end do
      end do
      opt = gbm_options(distribution=GBM_PAIRWISE, num_trees=10, interaction_depth=1, &
                        min_num_obs_in_node=2, shrinkage=0.08_dp, bag_fraction=0.65_dp, &
                        pairwise_metric=GBM_METRIC_NDCG)
      part = opt
      part%num_trees = 4
      call gbm_set_seed(7171)
      call gbm_fit(xr, yr, pair_full, opt, group=group)
      call gbm_set_seed(7171)
      call gbm_fit(xr, yr, pair_more, part, group=group)
      call gbm_continue(pair_more, xr, yr, 6, group=group)
      pr_full = gbm_predict(pair_full, xr)
      pr_more = gbm_predict(pair_more, xr)
      call check(maxval(abs(pr_full - pr_more)) < 1.0e-12_dp, 'pairwise continuation matches one-shot fit')
      call gbm_set_seed(7272)
      call gbm_cross_validate(xr, yr, cv_pair, 4, opt, group=group)
      call check(all(ieee_is_finite(cv_pair%error)) .and. all(ieee_is_finite(cv_pair%fitted)), &
                 'pairwise cross-validation is finite')
      do g = 1, ng
         call check(all(pack(cv_pair%fold_id, group == g) == cv_pair%fold_id(1 + (g - 1) * ni)), &
                    'pairwise group remains in one CV fold')
      end do

      do i = 1, ns
         xs(i, 1) = -1.2_dp + 2.4_dp * real(i - 1, dp) / real(ns - 1, dp)
         surv(i, 1) = 6.0_dp - 1.1_dp * xs(i, 1) + 0.03_dp * real(mod(i, 4), dp)
         surv(i, 2) = merge(1.0_dp, 0.0_dp, mod(i, 5) /= 0)
         strata(i) = 1
      end do
      opt = gbm_options(distribution=GBM_COXPH, num_trees=8, interaction_depth=1, &
                        min_num_obs_in_node=2, shrinkage=0.05_dp, bag_fraction=0.7_dp)
      part = opt
      part%num_trees = 3
      call gbm_set_seed(7373)
      call gbm_fit(xs, surv, cox_full, opt, strata=strata)
      call gbm_set_seed(7373)
      call gbm_fit(xs, surv, cox_more, part, strata=strata)
      call gbm_continue(cox_more, xs, surv, 5, strata=strata)
      ps_full = gbm_predict(cox_full, xs)
      ps_more = gbm_predict(cox_more, xs)
      call check(maxval(abs(ps_full - ps_more)) < 1.0e-12_dp, 'Cox continuation matches one-shot fit')
      call gbm_set_seed(7474)
      call gbm_cross_validate(xs, surv, cv_cox, 3, opt, strata=strata)
      call check(all(ieee_is_finite(cv_cox%error)) .and. all(ieee_is_finite(cv_cox%fitted)), &
                 'Cox cross-validation is finite')
   end subroutine test_pairwise_and_cox_expanded_apis

end program test_gbm3
