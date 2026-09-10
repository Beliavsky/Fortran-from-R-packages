program test_grf
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
   use grf
   implicit none

   integer, parameter :: n = 80
   integer, parameter :: p = 3
   real(dp) :: x(n,p)
   real(dp) :: y(n)
   real(dp) :: w(n)
   real(dp) :: z(n)
   real(dp) :: ym(n,2)
   real(dp) :: wm(n,2)
   integer :: cls(n)
   integer :: event(n)
   real(dp) :: time(n)
   type(grf_options) :: opt
   type(grf_options) :: opt_group
   type(grf_options) :: opt_tune
   type(grf_forest) :: rf
   type(grf_forest) :: grouped_rf
   type(grf_forest) :: tuned_rf
   type(grf_forest) :: grouped_cf
   type(grf_forest) :: grouped_ivf
   type(grf_forest) :: grouped_llf
   type(grf_forest) :: mrf
   type(grf_forest) :: pf
   type(grf_forest) :: qf
   type(grf_forest) :: cf
   type(grf_forest) :: ivf
   type(grf_forest) :: sf
   type(grf_forest) :: csf
   type(grf_forest) :: grouped_csf
   type(grf_forest) :: csf_probability
   type(grf_forest) :: lmf
   type(grf_forest) :: maf
   type(grf_forest) :: llf
   type(grf_forest) :: ll_split_forest
   type(grf_boosted_forest) :: bf
   type(grf_tree) :: tree
   type(grf_forest) :: pair(2)
   type(grf_forest) :: merged
   type(grf_ate_result) :: ate
   type(grf_linear_result) :: linear
   type(grf_rate_result) :: rate
   real(dp), allocatable :: pred(:)
   real(dp), allocatable :: pred2(:,:)
   real(dp), allocatable :: pred3(:,:,:)
   real(dp), allocatable :: pred_var(:)
   real(dp), allocatable :: kernel(:,:)
   real(dp), allocatable :: score(:)
   real(dp), allocatable :: importance(:)
   real(dp), allocatable :: sx(:,:)
   real(dp), allocatable :: sy(:)
   real(dp), allocatable :: stau(:)
   real(dp), allocatable :: sm(:)
   real(dp), allocatable :: se(:)
   real(dp), allocatable :: ssx(:,:)
   real(dp), allocatable :: ssy(:)
   real(dp), allocatable :: scate(:)
   real(dp), allocatable :: scatep(:)
   integer, allocatable :: sw(:)
   integer, allocatable :: sevent(:)
   integer, allocatable :: leaf(:)
   integer, allocatable :: freq(:,:)
   integer :: info
   integer :: i
   integer :: failures
   real(dp) :: rmse
   real(dp) :: selected_lambda
   real(dp) :: corr_num
   real(dp) :: corr_den
   real(dp) :: qv(3)
   real(dp) :: priorities(n)
   real(dp) :: a(n,1)

   failures = 0
   do i = 1, n
      x(i,1) = -1.0_dp + 2.0_dp * real(i - 1,dp) / real(n - 1,dp)
      x(i,2) = sin(0.37_dp * real(i,dp))
      x(i,3) = cos(0.19_dp * real(i,dp))
      y(i) = 2.0_dp * x(i,1) + 0.25_dp * x(i,2)
      if (mod(i,2) == 0) then
         w(i) = 1.0_dp
      else
         w(i) = 0.0_dp
      end if
      z(i) = w(i)
      cls(i) = merge(2,1,y(i) > 0.0_dp)
      time(i) = 0.6_dp + 0.8_dp * real(i,dp) / real(n,dp) + 0.2_dp * w(i)
      event(i) = merge(1,0,mod(i,4) /= 0)
      ym(i,1) = y(i)
      ym(i,2) = -0.5_dp * y(i) + x(i,3)
      wm(i,1) = w(i)
      wm(i,2) = merge(1.0_dp,0.0_dp,mod(i,3) == 0)
      priorities(i) = 1.0_dp + x(i,1)
      a(i,1) = x(i,1)
   end do
   x(7,3) = ieee_value(0.0_dp, ieee_quiet_nan)

   opt%num_trees = 36
   opt%sample_fraction = 0.8_dp
   opt%min_node_size = 4
   opt%honesty = .true.
   opt%honesty_fraction = 0.5_dp
   opt%seed = 1234

   call regression_forest(x, y, rf, info, opt)
   call check(info == 0, 'regression training', failures)
   call predict_regression_forest(rf, x, pred)
   rmse = sqrt(sum((pred - y)**2) / real(n,dp))
   call check(rmse < 0.75_dp, 'regression predictive signal', failures)

   opt_group = opt
   opt_group%num_trees = 5
   opt_group%ci_group_size = 2
   opt_group%sample_fraction = 0.5_dp
   call regression_forest(x, y, grouped_rf, info, opt_group)
   call check(info == 0, 'grouped-tree regression training', failures)
   call check(size(grouped_rf%trees) == 6, 'ci tree count rounds to full groups', failures)
   call check(all(grouped_rf%trees(1)%inbag .eqv. grouped_rf%trees(2)%inbag), &
              'ci group trees share half-sample at sample_fraction 0.5', failures)
   call predict_regression_forest(grouped_rf, x(1:5,:), pred, variance=pred_var)
   call check(size(pred_var) == 5 .and. all(.not. ieee_is_nan(pred_var)) .and. all(pred_var >= 0.0_dp), &
              'regression little-bag variance', failures)
   opt_tune = grf_options()
   opt_tune%num_trees = 12
   opt_tune%seed = 2468
   opt_tune%tune_parameters = .true.
   opt_tune%tune_num_trees = 6
   opt_tune%tune_num_reps = 6
   opt_tune%tune_num_draws = 20
   call regression_forest(x, y, tuned_rf, info, opt_tune)
   call check(info == 0 .and. size(tuned_rf%trees) == 12 .and. .not. tuned_rf%options%tune_parameters, &
              'regression automatic tuning', failures)
   call check(tuned_rf%options%sample_fraction > 0.0_dp .and. tuned_rf%options%sample_fraction <= 0.5_dp .and. &
              tuned_rf%options%mtry >= 0, 'tuned regression option ranges', failures)

   call get_forest_weights(rf, x(1:5,:), kernel)
   call check(maxval(abs(sum(kernel,dim=2) - 1.0_dp)) < 1.0e-10_dp, 'forest weights normalize', failures)
   call get_tree(rf, 1, tree, info)
   call check(info == 0 .and. tree%n_nodes >= 1, 'get_tree', failures)
   call get_leaf_node(tree, x(1:5,:), leaf)
   call check(all(leaf >= 1), 'get_leaf_node', failures)
   call split_frequencies(rf, freq)
   call check(sum(freq) > 0, 'split_frequencies', failures)
   call variable_importance(rf, importance)
   call check(abs(sum(importance) - 1.0_dp) < 1.0e-10_dp, 'variable_importance normalizes', failures)

   call multi_regression_forest(x, ym, mrf, info, opt)
   call check(info == 0, 'multi regression training', failures)
   call predict_multi_regression_forest(mrf, x(1:6,:), pred2)
   call check(all(shape(pred2) == [6,2]), 'multi regression prediction shape', failures)

   call probability_forest(x, cls, pf, info, opt)
   call check(info == 0, 'probability training', failures)
   call predict_probability_forest(pf, x(1:6,:), pred2)
   call check(maxval(abs(sum(pred2,dim=2) - 1.0_dp)) < 1.0e-10_dp, 'probabilities normalize', failures)

   qv = [0.2_dp, 0.5_dp, 0.8_dp]
   call quantile_forest(x, y, qv, qf, info, opt)
   call check(info == 0, 'quantile training', failures)
   call predict_quantile_forest(qf, x(1:6,:), pred2)
   call check(all(pred2(:,1) <= pred2(:,2)) .and. all(pred2(:,2) <= pred2(:,3)), 'quantiles ordered', failures)

   y = 0.5_dp * x(:,1) + (1.0_dp + 0.8_dp * x(:,1)) * w
   call causal_forest(x, y, w, cf, info, opt)
   call check(info == 0, 'causal training', failures)
   call predict_causal_forest(cf, x, pred)
   call check(count(.not. ieee_is_nan(pred)) > n / 2, 'causal predictions finite', failures)
   corr_num = sum((pred - sum(pred)/real(n,dp)) * (priorities - sum(priorities)/real(n,dp)), &
                  mask=.not. ieee_is_nan(pred))
   corr_den = sqrt(sum((pred - sum(pred)/real(n,dp))**2, mask=.not. ieee_is_nan(pred)) * &
                   sum((priorities - sum(priorities)/real(n,dp))**2))
   if (corr_den > 0.0_dp) call check(corr_num / corr_den > 0.05_dp, 'causal heterogeneity signal', failures)
   opt_group%num_trees = 12
   call causal_forest(x, y, w, grouped_cf, info, opt_group)
   call check(info == 0, 'grouped causal training', failures)
   call predict_causal_forest(grouped_cf, x(1:6,:), pred, variance=pred_var)
   call check(count(.not. ieee_is_nan(pred_var)) >= 4 .and. all((pred_var >= 0.0_dp) .or. ieee_is_nan(pred_var)), &
              'causal little-bag variance', failures)
   call get_scores_causal_forest(cf, score)
   call check(size(score) == n, 'causal score length', failures)
   call average_treatment_effect(cf, ate, info)
   call check(info == 0 .and. size(ate%estimate) == 1, 'average_treatment_effect', failures)
   call test_calibration(cf, linear, info)
   call check(info == 0, 'test_calibration', failures)
   call best_linear_projection(cf, a, linear, info)
   call check(info == 0 .and. size(linear%coefficient) == 2, 'best_linear_projection', failures)
   call rank_average_treatment_effect_fit(1.0_dp + 0.8_dp * x(:,1), priorities, rate, info, bootstrap_reps=20, seed=9)
   call check(info == 0 .and. rate%estimate > 0.0_dp, 'RATE fit', failures)
   call rank_average_treatment_effect(cf, priorities, rate, info, bootstrap_reps=10, seed=9)
   call check(info == 0, 'RATE forest', failures)

   call instrumental_forest(x, y, w, z, ivf, info, opt)
   call check(info == 0, 'instrumental training', failures)
   call check(allocated(ivf%compliance_hat) .and. count(abs(ivf%compliance_hat) > 1.0e-8_dp) > n / 2, &
              'instrumental compliance nuisance retained', failures)
   call predict_instrumental_forest(ivf, x(1:6,:), pred)
   call check(size(pred) == 6, 'instrumental prediction shape', failures)
   call get_scores_instrumental_forest(ivf, score)
   call check(size(score) == n .and. count(.not. ieee_is_nan(score)) > n / 2, &
              'instrumental doubly robust scores', failures)
   call instrumental_forest(x, y, w, z, grouped_ivf, info, opt_group)
   call check(info == 0, 'grouped instrumental training', failures)
   call predict_instrumental_forest(grouped_ivf, x(1:6,:), pred, variance=pred_var)
   call check(count(.not. ieee_is_nan(pred_var)) >= 4 .and. all((pred_var >= 0.0_dp) .or. ieee_is_nan(pred_var)), &
              'instrumental little-bag variance', failures)

   call survival_forest(x, time, event, sf, info, opt)
   call check(info == 0, 'survival training', failures)
   call predict_survival_forest(sf, x(1:4,:), pred2)
   call check(size(pred2,1) == 4 .and. all(pred2 >= 0.0_dp) .and. all(pred2 <= 1.0_dp), 'survival predictions', failures)

   call causal_survival_forest(x, time, w, event, 1.0_dp, csf, info, opt)
   call check(info == 0, 'causal survival training', failures)
   call check(allocated(csf%causal_survival_numerator) .and. allocated(csf%causal_survival_denominator), &
              'causal survival doubly robust moments retained', failures)
   call check(allocated(csf%censor_survival_at_y) .and. minval(csf%censor_survival_at_y) > 0.0_dp, &
              'causal survival censoring nuisance retained', failures)
   call check(allocated(csf%treatment_variance_hat) .and. all(csf%treatment_variance_hat >= 0.0_dp), &
              'causal survival treatment variance nuisance retained', failures)
   call predict_causal_survival_forest(csf, x(1:5,:), pred)
   call check(size(pred) == 5, 'causal survival prediction shape', failures)
   call get_scores_causal_survival_forest(csf, score)
   call check(size(score) == n .and. count(.not. ieee_is_nan(score)) > n / 2, &
              'causal survival doubly robust scores', failures)
   call causal_survival_forest(x, time, w, event, 1.0_dp, grouped_csf, info, opt_group)
   call check(info == 0, 'grouped causal survival training', failures)
   call predict_causal_survival_forest(grouped_csf, x(1:5,:), pred, variance=pred_var)
   call check(count(.not. ieee_is_nan(pred_var)) >= 3 .and. &
              all((pred_var >= 0.0_dp) .or. ieee_is_nan(pred_var)), &
              'causal survival little-bag variance', failures)
   call causal_survival_forest(x, time, w, event, 1.0_dp, csf_probability, info, opt, &
                               survival_probability=.true.)
   call check(info == 0 .and. csf_probability%survival_probability_target, &
              'causal survival probability target training', failures)
   call predict_causal_survival_forest(csf_probability, x(1:5,:), pred)
   call check(size(pred) == 5 .and. count(.not. ieee_is_nan(pred)) >= 3, &
              'causal survival probability predictions', failures)

   call lm_forest(x, ym, wm, lmf, info, opt)
   call check(info == 0, 'lm forest training', failures)
   call predict_lm_forest(lmf, x(1:4,:), pred3)
   call check(all(shape(pred3) == [4,2,2]), 'lm forest prediction shape', failures)

   call multi_arm_causal_forest(x, ym(:,1), wm, maf, info, opt)
   call check(info == 0, 'multi-arm forest training', failures)
   call predict_multi_arm_causal_forest(maf, x(1:4,:), pred2)
   call check(all(shape(pred2) == [4,2]), 'multi-arm prediction shape', failures)

   call ll_regression_forest(x, ym(:,1), llf, info, opt)
   call check(info == 0, 'll regression training', failures)
   call predict_ll_regression_forest(llf, x(1:4,:), pred, ridge=1.0e-4_dp)
   call check(size(pred) == 4 .and. count(.not. ieee_is_nan(pred)) >= 3, 'll regression prediction', failures)
   call predict_ll_regression_forest(llf, x(1:4,:), pred, selected_ridge=selected_lambda)
   call check(selected_lambda >= 0.0_dp .and. selected_lambda <= 10.0_dp .and. &
              count(.not. ieee_is_nan(pred)) >= 3, 'll regression OOB ridge tuning', failures)
   call ll_regression_forest(x, ym(:,1), ll_split_forest, info, opt, enable_ll_split=.true., &
                             split_ridge=0.1_dp, split_cutoff=8)
   call check(info == 0 .and. ll_split_forest%ll_split_enabled .and. &
              size(ll_split_forest%ll_split_variables) == p, 'local-linear residual split training', failures)
   call predict_ll_regression_forest(ll_split_forest, x(1:4,:), pred, ridge=0.1_dp, weight_penalty=.true.)
   call check(count(.not. ieee_is_nan(pred)) >= 3, 'local-linear residual split prediction', failures)
   call ll_regression_forest(x, ym(:,1), grouped_llf, info, opt_group)
   call check(info == 0, 'grouped ll regression training', failures)
   call predict_ll_regression_forest(grouped_llf, x(1:4,:), pred, ridge=1.0e-4_dp, variance=pred_var)
   call check(count(.not. ieee_is_nan(pred_var)) >= 3 .and. all((pred_var >= 0.0_dp) .or. ieee_is_nan(pred_var)), &
              'local-linear little-bag variance', failures)

   call boosted_regression_forest(x, ym(:,1), bf, info, num_stages=3, learning_rate=0.2_dp, options=opt)
   call check(info == 0, 'boosted regression training', failures)
   call check(bf%n_stages == 3 .and. size(bf%stage_error) == 3 .and. &
              allocated(bf%training_prediction), 'boosted fixed-stage diagnostics', failures)
   call predict_boosted_regression_forest(bf, x(1:4,:), pred)
   call check(size(pred) == 4, 'boosted prediction shape', failures)
   call boosted_regression_forest(x, ym(:,1), bf, info, options=opt, error_reduction=0.99_dp, &
                                  max_stages=3, tune_trees=8)
   call check(info == 0 .and. bf%n_stages >= 1 .and. bf%n_stages <= 3, &
              'boosted automatic stage selection', failures)
   call check(size(bf%stage_error) == bf%n_stages .and. &
              count(.not. ieee_is_nan(bf%training_prediction)) > n / 2, &
              'boosted OOB debiased error diagnostics', failures)

   pair(1) = rf
   pair(2) = rf
   call merge_forests(pair, merged, info)
   call check(info == 0 .and. size(merged%trees) == 2 * size(rf%trees), 'merge_forests', failures)

   call generate_causal_data(30, 3, sx, sy, sw, stau, sm, se, info, seed=11)
   call check(info == 0 .and. all(shape(sx) == [30,3]) .and. size(sy) == 30, 'generate_causal_data', failures)
   call generate_causal_survival_data(30, 2, ssx, ssy, sw, sevent, scate, scatep, info, seed=12)
   call check(info == 0 .and. all(shape(ssx) == [30,2]) .and. size(ssy) == 30, 'generate_causal_survival_data', failures)

   if (failures /= 0) then
      error stop 'grf tests failed'
   end if
   print '(a)', 'All grf deterministic tests passed.'

contains

   subroutine check(condition, label, failures)
      logical, intent(in) :: condition !! Boolean assertion result to record.
      character(len=*), intent(in) :: label !! Short deterministic test label printed on failure.
      integer, intent(inout) :: failures !! Running number of failed assertions.

      if (.not. condition) then
         failures = failures + 1
         print '(a,a)', 'FAIL: ', trim(label)
      end if
   end subroutine check

end program test_grf
