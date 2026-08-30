! SPDX-License-Identifier: GPL-3.0-only
! See NOTICE.md for ranger upstream provenance.
program test_ranger
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp, i64
   use ranger, only : ranger_options
   use ranger, only : ranger_classification_forest, ranger_probability_forest
   use ranger, only : ranger_regression_forest, ranger_survival_forest
   use ranger, only : fit_ranger_classification, predict_ranger_classification
   use ranger, only : fit_ranger_probability, predict_ranger_probability
   use ranger, only : fit_ranger_regression, predict_ranger_regression, predict_ranger_quantiles, predict_ranger_quantiles_oob
   use ranger, only : fit_ranger_survival, predict_ranger_survival
   use ranger, only : RANGER_IMPORTANCE_PERMUTATION, RANGER_SPLIT_HELLINGER, RANGER_SPLIT_EXTRATREES
   use ranger, only : RANGER_SPLIT_BETA, RANGER_SPLIT_POISSON
   use ranger, only : RANGER_UNORDERED_PARTITION, RANGER_UNORDERED_ORDER, RANGER_IMPORTANCE_IMPURITY_CORRECTED
   use ranger, only : importance_pvalues_janitza, importance_pvalues_altmann, infinitesimal_jackknife
   use ranger, only : hierarchical_shrink_regression, hierarchical_shrink_probability
   use ranger, only : case_specific_weights, tree_sizes, variable_usage
   use ranger_rng, only : ranger_rng_state, sample_indices, shuffle_int
   use ranger_maxstat, only : average_ranks, logrank_scores, maxstat_best
   use ranger_maxstat, only : maxstat_pvalue_lau92, maxstat_pvalue_lau94, maxstat_pvalue_unadjusted
   use ranger_maxstat, only : adjust_pvalues_bh
   use ranger_forest_common, only : draw_inbag, draw_classwise_inbag, validate_options, tree_options_for_index
   use ranger_tree_common, only : candidate_variables
   use ranger_factor_order, only : order_factors_survival
   use ranger_tree_survival, only : concordance_index, concordance_casewise
   use ranger_types, only : RANGER_TERMINAL
   implicit none

   call test_rng_parity()
   call test_sampling_parity()
   call test_candidate_parity()
   call test_option_parity()
   call test_oob_error_option()
   call test_oob_no_coverage_parity()
   call test_maxstat_parity()
   call test_classification()
   call test_classification_vote_ties()
   call test_categorical_missing_holdout()
   call test_probability()
   call test_regression()
   call test_special_regression_rules()
   call test_bootstrap_duplicate_semantics()
   call test_corrected_impurity()
   call test_ordered_factor_preprocessing()
   call test_survival_factor_ordering()
   call test_survival_time_interest_parity()
   call test_oob_quantiles()
   call test_survival()
   call test_survival_concordance_parity()
   call test_utilities()
   call test_ij_small_sample()
   call test_ij_calibration()
   print '(a)', 'all ranger tests passed'

contains

   subroutine test_rng_parity()
      type(ranger_rng_state) :: rng
      real(dp) :: expected(5), weights(5)
      integer :: i, sample5(5), sample3(3)

      expected = [0.78682095486780201_dp, 0.25048034068802871_dp, 0.71067122897865553_dp, &
         0.94666780096097036_dp, 0.019271058195813772_dp]
      call rng%seed(5489_i64)
      do i = 1, 5
         call require(abs(rng%uniform() - expected(i)) < 2.0e-16_dp, 'mt19937_64 uniform stream')
      end do
      call rng%seed(5489_i64)
      call sample_indices(rng, 100, 5, .false., sample5)
      call require(all(sample5 == [79, 26, 72, 95, 2]), 'ranger simple sampling stream')
      call rng%seed(5489_i64)
      call sample_indices(rng, 10, 5, .false., sample5)
      call require(all(sample5 == [8, 4, 1, 10, 5]), 'ranger Fisher-Yates sampling stream')
      weights = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
      call rng%seed(5489_i64)
      call sample_indices(rng, 5, 3, .false., sample3, weights=weights)
      call require(all(sample3 == [5, 3, 1]), 'ranger weighted sampling stream')
   end subroutine test_rng_parity

   subroutine test_sampling_parity()
      type(ranger_rng_state) :: rng
      type(ranger_options) :: options
      integer :: perm(10), inbag(10), y(10), i, status
      real(dp) :: case_weight(10), u

      call rng%seed(42_i64)
      perm = [(i, i = 1, 10)]
      call shuffle_int(rng, perm)
      call require(all(perm == [7, 10, 2, 5, 6, 4, 1, 8, 9, 3]), 'libstdc++ std::shuffle stream')

      options%replace = .false.
      options%sample_fraction = 0.5_dp
      case_weight = 1.0_dp
      call rng%seed(42_i64)
      call draw_inbag(rng, options, case_weight, .false., inbag, status)
      call require(status == 0, 'unweighted no-replacement sampling status')
      call require(all(inbag == [0, 1, 0, 0, 1, 1, 1, 0, 0, 1]), &
         'shuffleAndSplit in-bag selection')
      u = rng%uniform()
      call require(abs(u - 0.75515553295453897_dp) < 2.0e-16_dp, &
         'shuffleAndSplit RNG pass-by-value semantics')

      y = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2]
      options%sample_fraction = 1.0_dp
      call rng%seed(42_i64)
      call draw_classwise_inbag(rng, y, 2, options, [0.2_dp, 0.2_dp], inbag, status)
      call require(status == 0, 'classwise no-replacement sampling status')
      u = rng%uniform()
      call require(abs(u - 0.75515553295453897_dp) < 2.0e-16_dp, &
         'classwise shuffle RNG pass-by-value semantics')
   end subroutine test_sampling_parity

   subroutine test_candidate_parity()
      type(ranger_rng_state) :: rng
      type(ranger_options) :: options
      integer, allocatable :: vars(:)
      integer :: nvars

      options%importance_mode = RANGER_IMPORTANCE_IMPURITY_CORRECTED
      options%mtry = 2
      allocate(options%always_split_variables(1), options%split_select_weights(4))
      options%always_split_variables = [4]
      options%split_select_weights = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp]
      call rng%seed(3_i64)
      call candidate_variables(rng, 8, 2, options, vars, nvars)
      call require(nvars == 4, 'corrected-importance candidate count')
      call require(all(vars == [5, 2, 4, 5]), 'corrected-importance weighted shadow candidate semantics')
   end subroutine test_candidate_parity

   subroutine test_option_parity()
      type(ranger_options) :: options, tree_options
      integer :: rc
      character(len=200) :: message
      real(dp) :: x(5,1)
      integer :: y(5)
      type(ranger_classification_forest) :: forest
      type(ranger_probability_forest) :: probability_forest

      options%num_trees = 2
      options%mtry = 1
      allocate(options%split_select_weights_by_tree(2, 2))
      options%split_select_weights_by_tree = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
      call validate_options(options, 2, rc, message)
      call require(rc == 0, 'per-tree split-select weights validate')
      call tree_options_for_index(options, 1, tree_options)
      call require(allocated(tree_options%split_select_weights), 'per-tree split weights materialized')
      call require(maxval(abs(tree_options%split_select_weights - options%split_select_weights_by_tree(1, :))) < 1.0e-15_dp, &
         'per-tree split-select first row')
      call tree_options_for_index(options, 2, tree_options)
      call require(maxval(abs(tree_options%split_select_weights - options%split_select_weights_by_tree(2, :))) < 1.0e-15_dp, &
         'per-tree split-select second row')

      if (allocated(options%split_select_weights_by_tree)) deallocate(options%split_select_weights_by_tree)
      allocate(options%regularization_factor(1))
      options%regularization_factor = [0.8_dp]
      call validate_options(options, 2, rc, message)
      call require(rc == 0, 'scalar regularization factor accepted')

      x(:,1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
      y = [1, 1, 1, 2, 2]
      options = ranger_options()
      options%num_trees = 1
      options%mtry = 1
      options%replace = .false.
      call fit_ranger_classification(x, y, forest, options=options, &
         sample_fraction_class=[0.61_dp, 0.39_dp], status=rc, message=message)
      call require(rc /= 0, 'class-wise no-replacement raw sample-fraction validation')

      ! The R interface accepts nonnegative class-specific limits.  A scalar
      ! zero is ranger's default sentinel, but zeros in a vector are retained.
      options = ranger_options()
      options%num_trees = 1
      options%mtry = 1
      options%replace = .false.
      call fit_ranger_classification(x, y, forest, options=options, &
         min_node_size_class=[0, 0], min_bucket_class=[0, 0], status=rc, message=message)
      call require(rc == 0, 'zero class-specific node limits accepted for classification')
      call fit_ranger_probability(x, y, probability_forest, options=options, &
         min_node_size_class=[0, 0], min_bucket_class=[0, 0], status=rc, message=message)
      call require(rc == 0, 'zero class-specific node limits accepted for probability')
   end subroutine test_option_parity

   subroutine test_oob_error_option()
      real(dp) :: x(12, 1), yreg(12), time(12)
      integer :: yclass(12), event(12), i, status
      type(ranger_options) :: options
      type(ranger_classification_forest) :: class_forest
      type(ranger_probability_forest) :: probability_forest
      type(ranger_regression_forest) :: regression_forest
      type(ranger_survival_forest) :: survival_forest

      do i = 1, 12
         x(i, 1) = real(i, dp)
         yreg(i) = 0.5_dp * real(i, dp)
         time(i) = real(13 - i, dp)
         yclass(i) = merge(1, 2, i <= 6)
      end do
      event = 1
      options%num_trees = 5
      options%mtry = 1
      options%min_node_size = 1
      options%min_bucket = 1
      options%importance_mode = RANGER_IMPORTANCE_PERMUTATION
      options%oob_error = .false.
      options%seed = 24601_i64

      call fit_ranger_classification(x, yclass, class_forest, options=options, status=status)
      call require(status == 0, 'classification oob.error false fit')
      call require(all(class_forest%oob_count == 0), 'classification oob.error false skips OOB prediction')
      call require(all(ieee_is_nan(class_forest%prediction_error)), 'classification oob.error false NaN error')

      call fit_ranger_probability(x, yclass, probability_forest, options=options, status=status)
      call require(status == 0, 'probability oob.error false fit')
      call require(all(probability_forest%oob_count == 0), 'probability oob.error false skips OOB prediction')
      call require(all(ieee_is_nan(probability_forest%prediction_error)), 'probability oob.error false NaN error')

      call fit_ranger_regression(x, yreg, regression_forest, options=options, status=status)
      call require(status == 0, 'regression oob.error false fit')
      call require(all(regression_forest%oob_count == 0), 'regression oob.error false skips OOB prediction')
      call require(all(ieee_is_nan(regression_forest%prediction_error)), 'regression oob.error false NaN error')

      call fit_ranger_survival(x, time, event, survival_forest, options=options, status=status)
      call require(status == 0, 'survival oob.error false fit')
      call require(all(survival_forest%oob_count == 0), 'survival oob.error false skips OOB prediction')
      call require(all(ieee_is_nan(survival_forest%prediction_error)), 'survival oob.error false NaN error')
   end subroutine test_oob_error_option

   subroutine test_oob_no_coverage_parity()
      real(dp) :: x(12, 1), yreg(12), time(12)
      integer :: yclass(12), event(12), i, status
      type(ranger_options) :: options
      type(ranger_classification_forest) :: class_forest
      type(ranger_probability_forest) :: probability_forest
      type(ranger_regression_forest) :: regression_forest
      type(ranger_survival_forest) :: survival_forest

      do i = 1, 12
         x(i, 1) = real(i, dp)
         yreg(i) = real(i, dp) / 3.0_dp
         time(i) = real(i, dp)
         yclass(i) = merge(1, 2, i <= 6)
      end do
      event = 1
      options%num_trees = 3
      options%mtry = 1
      options%min_node_size = 1
      options%min_bucket = 1
      options%replace = .false.
      options%sample_fraction = 1.0_dp
      options%seed = 271828_i64

      call fit_ranger_classification(x, yclass, class_forest, options=options, status=status)
      call require(status == 0, 'classification zero-OOB fit')
      call require(all(class_forest%oob_count == 0), 'classification zero-OOB counts')
      call require(all(ieee_is_nan(class_forest%prediction_error)), 'classification zero-OOB error is NaN')

      call fit_ranger_probability(x, yclass, probability_forest, options=options, status=status)
      call require(status == 0, 'probability zero-OOB fit')
      call require(all(probability_forest%oob_count == 0), 'probability zero-OOB counts')
      call require(all(ieee_is_nan(probability_forest%prediction_error)), 'probability zero-OOB error is NaN')
      call require(all(ieee_is_nan(probability_forest%oob_probability)), 'probability zero-OOB predictions are NaN')

      call fit_ranger_regression(x, yreg, regression_forest, options=options, status=status)
      call require(status == 0, 'regression zero-OOB fit')
      call require(all(regression_forest%oob_count == 0), 'regression zero-OOB counts')
      call require(all(ieee_is_nan(regression_forest%prediction_error)), 'regression zero-OOB error is NaN')
      call require(all(ieee_is_nan(regression_forest%oob_prediction)), 'regression zero-OOB predictions are NaN')

      call fit_ranger_survival(x, time, event, survival_forest, options=options, status=status)
      call require(status == 0, 'survival zero-OOB fit')
      call require(all(survival_forest%oob_count == 0), 'survival zero-OOB counts')
      call require(all(ieee_is_nan(survival_forest%prediction_error)), 'survival zero-OOB error is NaN')
      call require(all(abs(survival_forest%oob_chf) <= 0.0_dp), 'survival zero-OOB CHF remains zero')
   end subroutine test_oob_no_coverage_parity

   subroutine test_maxstat_parity()
      real(dp) :: values(4), ranks(4), times(4), scores(4), p(3), adjusted(3)
      real(dp) :: xcut(4), best_stat, best_split
      integer, allocatable :: num_left(:)
      integer :: status4(4), cuts(4)
      logical :: found

      values = [1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp]
      call average_ranks(values, ranks)
      call require(maxval(abs(ranks - [1.0_dp, 2.5_dp, 2.5_dp, 4.0_dp])) < 1.0e-14_dp, 'average ranks')
      times = [1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
      status4 = [1, 0, 1, 1]
      call logrank_scores(times, status4, scores)
      call require(maxval(abs(scores - [2.0_dp / 3.0_dp, -1.0_dp / 3.0_dp, 1.0_dp / 6.0_dp, &
         -5.0_dp / 6.0_dp])) < 1.0e-14_dp, 'logrank scores')
      call require(abs(maxstat_pvalue_lau92(2.0_dp, 0.1_dp, 0.9_dp) - 0.463872768757117_dp) < 1.0e-14_dp, &
         'Lau92 p-value')
      cuts = [2, 5, 8, 10]
      call require(abs(maxstat_pvalue_lau94(2.0_dp, 0.1_dp, 0.9_dp, 12, cuts) - &
         0.14996182821806997_dp) < 1.0e-14_dp, 'Lau94 p-value')
      p = [0.04_dp, 0.01_dp, 0.03_dp]
      call adjust_pvalues_bh(p, adjusted)
      call require(maxval(abs(adjusted - [0.04_dp, 0.03_dp, 0.04_dp])) < 1.0e-14_dp, 'BH adjustment')

      xcut = [1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
      call maxstat_best(ranks, xcut, 0.1_dp, 0.9_dp, best_stat, best_split, num_left, found)
      call require(found, 'maxstat cutpoint found')
      call require(all(num_left == [2, 3, 4]), 'maxstat cumulative cutpoint counts')
      call require(abs(maxstat_pvalue_unadjusted(2.0_dp) - erfc(sqrt(2.0_dp))) < 1.0e-15_dp, &
         'unadjusted maxstat p-value')
   end subroutine test_maxstat_parity

   subroutine test_classification()
      real(dp) :: x(24, 2), votes(24, 2)
      integer :: y(24), prediction(24), i, status
      type(ranger_options) :: options
      type(ranger_classification_forest) :: forest

      do i = 1, 12
         x(i, 1) = -real(13 - i, dp)
         x(i, 2) = real(mod(i, 3) + 1, dp)
         y(i) = 1
      end do
      do i = 13, 24
         x(i, 1) = real(i - 12, dp)
         x(i, 2) = real(mod(i, 3) + 1, dp)
         y(i) = 2
      end do
      options%num_trees = 31
      options%mtry = 2
      options%importance_mode = RANGER_IMPORTANCE_PERMUTATION
      options%keep_inbag = .true.
      options%seed = 12345
      call fit_ranger_classification(x, y, forest, options=options, status=status)
      call require(status == 0, 'classification fit status')
      call predict_ranger_classification(forest, x, prediction, votes)
      call require(count(prediction == y) >= 23, 'classification training predictions')
      call require(all(abs(sum(votes, dim=2) - 1.0_dp) < 1.0e-12_dp), 'classification vote normalization')
      call require(allocated(forest%inbag), 'classification inbag retained')
      call require(count(forest%oob_count > 0) >= 20, 'classification OOB coverage')
      call require(forest%variable_importance(1) >= forest%variable_importance(2) - 0.05_dp, &
         'classification permutation importance')

      options%split_rule = RANGER_SPLIT_HELLINGER
      options%num_trees = 3
      options%importance_mode = 0
      call fit_ranger_classification(x, y, forest, options=options, status=status)
      call require(status == 0, 'Hellinger fit status')
      call predict_ranger_classification(forest, x, prediction)
      call require(count(prediction == y) >= 22, 'Hellinger predictions')
   end subroutine test_classification

   subroutine test_classification_vote_ties()
      type(ranger_classification_forest) :: forest
      real(dp) :: x(2, 1), votes(2, 2)
      integer :: prediction(2), t

      forest%nclass = 2
      forest%nvar = 1
      forest%nobs = 0
      allocate(forest%ncat(1), forest%category_map(1, 1), forest%trees(2))
      forest%ncat = 1
      forest%category_map = 1
      do t = 1, 2
         forest%trees(t)%n_nodes = 1
         allocate(forest%trees(t)%status(1), forest%trees(t)%node_class(1))
         forest%trees(t)%status = RANGER_TERMINAL
         forest%trees(t)%node_class = t
      end do
      x = 0.0_dp
      call predict_ranger_classification(forest, x, prediction, vote_fraction=votes, seed=42_i64)
      call require(all(abs(votes - 0.5_dp) < 1.0e-15_dp), 'classification tied vote fractions')
      call require(all(prediction == prediction(1)), 'classification tie RNG is passed by value')
      call require(prediction(1) == 2, 'classification tie uses mt19937_64 uniform choice')
   end subroutine test_classification_vote_ties


   subroutine test_categorical_missing_holdout()
      real(dp) :: x(20, 2), case_weight(20), nan_value
      integer :: y(20), ncat(2), prediction(20), i, status
      type(ranger_options) :: options
      type(ranger_classification_forest) :: forest

      do i = 1, 10
         x(i, 1) = 1.0_dp
         x(i, 2) = real(mod(i, 4), dp)
         y(i) = 1
      end do
      do i = 11, 20
         x(i, 1) = 2.0_dp
         x(i, 2) = real(mod(i, 4), dp)
         y(i) = 2
      end do
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      x(2, 1) = nan_value
      ncat = [2, 1]
      options%num_trees = 7
      options%mtry = 2
      options%replace = .false.
      options%sample_fraction = 1.0_dp
      options%respect_unordered_factors = RANGER_UNORDERED_PARTITION
      options%seed = 4321
      call fit_ranger_classification(x, y, forest, ncat=ncat, options=options, status=status)
      call require(status == 0, 'categorical/missing fit status')
      call predict_ranger_classification(forest, x, prediction)
      call require(count(prediction == y) >= 19, 'categorical/missing predictions')

      case_weight = 1.0_dp
      do i = 1, 20, 3
         case_weight(i) = 0.0_dp
      end do
      options%holdout = .true.
      options%replace = .true.
      options%sample_fraction = 1.0_dp
      call fit_ranger_classification(x, y, forest, ncat=ncat, options=options, case_weights=case_weight, status=status)
      call require(status == 0, 'holdout fit status')
      do i = 1, 20
         if (case_weight(i) <= 0.0_dp) then
            call require(forest%oob_count(i) > 0, 'holdout observations evaluated')
         else
            call require(forest%oob_count(i) == 0, 'training observations excluded from holdout error')
         end if
      end do
   end subroutine test_categorical_missing_holdout

   subroutine test_probability()
      real(dp) :: x(20, 2), probability(20, 2), before, after
      integer :: y(20), prediction(20), i, status
      type(ranger_options) :: options
      type(ranger_probability_forest) :: forest

      do i = 1, 10
         x(i, :) = [-real(11 - i, dp), real(mod(i, 2) + 1, dp)]
         y(i) = 1
      end do
      do i = 11, 20
         x(i, :) = [real(i - 10, dp), real(mod(i, 2) + 1, dp)]
         y(i) = 2
      end do
      options%num_trees = 25
      options%mtry = 2
      options%min_node_size = 1
      options%seed = 9876
      call fit_ranger_probability(x, y, forest, options=options, status=status)
      call require(status == 0, 'probability fit status')
      call predict_ranger_probability(forest, x, probability, prediction)
      call require(count(prediction == y) >= 19, 'probability predictions')
      call require(all(abs(sum(probability, dim=2) - 1.0_dp) < 1.0e-12_dp), 'probability normalization')
      before = forest%trees(1)%class_prob(1, forest%trees(1)%left(1))
      call hierarchical_shrink_probability(forest, 5.0_dp)
      after = forest%trees(1)%class_prob(1, forest%trees(1)%left(1))
      call require(ieee_is_finite(before) .and. ieee_is_finite(after), 'probability hierarchical shrinkage')
   end subroutine test_probability

   subroutine test_regression()
      real(dp) :: x(40, 2), y(40), prediction(40), mse, q(4, 3), probs(3)
      real(dp) :: root_before, leaf_after
      integer :: i, status
      type(ranger_options) :: options
      type(ranger_regression_forest) :: forest

      do i = 1, 40
         x(i, 1) = -2.0_dp + 4.0_dp * real(i - 1, dp) / 39.0_dp
         x(i, 2) = sin(real(i, dp))
         y(i) = 2.0_dp * x(i, 1) + 0.2_dp * x(i, 2)
      end do
      options%num_trees = 41
      options%mtry = 2
      options%min_node_size = 2
      options%quantreg = .true.
      options%keep_inbag = .false.
      options%importance_mode = RANGER_IMPORTANCE_PERMUTATION
      options%seed = 777
      call fit_ranger_regression(x, y, forest, options=options, status=status)
      call require(status == 0, 'regression fit status')
      call predict_ranger_regression(forest, x, prediction)
      mse = sum((prediction - y) ** 2) / real(size(y), dp)
      call require(mse < 0.35_dp, 'regression training MSE')
      probs = [0.1_dp, 0.5_dp, 0.9_dp]
      call predict_ranger_quantiles(forest, x([5, 15, 25, 35], :), probs, q)
      call require(all(q(:, 1) <= q(:, 2)) .and. all(q(:, 2) <= q(:, 3)), 'quantile ordering')
      root_before = forest%trees(1)%node_mean(1)
      call hierarchical_shrink_regression(forest, 10.0_dp)
      leaf_after = forest%trees(1)%node_mean(first_terminal(forest%trees(1)%status, forest%trees(1)%n_nodes))
      call require(ieee_is_finite(root_before) .and. ieee_is_finite(leaf_after), 'regression hierarchical shrinkage')
   end subroutine test_regression

   subroutine test_special_regression_rules()
      real(dp) :: x(30, 1), y(30), pred(30)
      integer :: i, status
      type(ranger_options) :: options
      type(ranger_regression_forest) :: forest

      do i = 1, 30
         x(i, 1) = real(i, dp) / 31.0_dp
         y(i) = 0.1_dp + 0.8_dp * x(i, 1)
      end do
      options%num_trees = 3
      options%mtry = 1
      options%min_node_size = 3
      options%replace = .false.
      options%sample_fraction = 1.0_dp
      y = [(real(i - 1, dp) / 29.0_dp, i = 1, 30)]
      options%split_rule = RANGER_SPLIT_BETA
      call fit_ranger_regression(x, y, forest, options=options, status=status)
      call require(status == 0, 'beta split rule including boundary responses')
      call predict_ranger_regression(forest, x, pred)
      call require(all(ieee_is_finite(pred)), 'beta predictions finite')

      y = real([(mod(i, 5), i = 1, 30)], dp)
      options%split_rule = RANGER_SPLIT_POISSON
      call fit_ranger_regression(x, y, forest, options=options, status=status)
      call require(status == 0, 'Poisson split rule')
      y = 0.0_dp
      call fit_ranger_regression(x, y, forest, options=options, status=status)
      call require(status /= 0, 'Poisson rejects zero-total response')

      y = 2.0_dp * x(:, 1)
      options%split_rule = RANGER_SPLIT_EXTRATREES
      options%num_random_splits = 8
      call fit_ranger_regression(x, y, forest, options=options, status=status)
      call require(status == 0, 'ExtraTrees split rule')
   end subroutine test_special_regression_rules

   subroutine test_bootstrap_duplicate_semantics()
      real(dp) :: x(4, 1)
      integer :: y(4), preset(4, 1), status
      type(ranger_options) :: options
      type(ranger_classification_forest) :: forest

      x(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
      y = [1, 2, 1, 2]
      preset(:, 1) = [2, 2, 0, 0]
      options%num_trees = 1
      options%mtry = 1
      options%min_node_size = 3
      options%min_bucket = 1
      options%seed = 5_i64
      call fit_ranger_classification(x, y, forest, options=options, preset_inbag=preset, status=status)
      call require(status == 0, 'duplicate-bootstrap fit status')
      call require(forest%trees(1)%node_n(1) == 4, 'bootstrap multiplicities count toward node size')
      call require(forest%trees(1)%n_nodes == 3, 'duplicate bootstrap permits ranger root split')
   end subroutine test_bootstrap_duplicate_semantics

   subroutine test_corrected_impurity()
      real(dp) :: x(40, 3)
      integer :: y(40), i, status
      type(ranger_options) :: options
      type(ranger_classification_forest) :: forest

      do i = 1, 40
         x(i, 1) = real(i, dp)
         x(i, 2) = real(mod(7 * i, 11), dp)
         x(i, 3) = sin(real(i, dp))
         y(i) = merge(1, 2, i <= 20)
      end do
      options%num_trees = 61
      options%mtry = 3
      options%min_node_size = 2
      options%importance_mode = RANGER_IMPORTANCE_IMPURITY_CORRECTED
      options%seed = 1234_i64
      call fit_ranger_classification(x, y, forest, options=options, status=status)
      call require(status == 0, 'corrected impurity fit')
      call require(all(ieee_is_finite(forest%variable_importance)), 'corrected impurity finite')
      call require(any(forest%variable_importance < 0.0_dp), 'corrected impurity shadow subtraction')
   end subroutine test_corrected_impurity

   subroutine test_ordered_factor_preprocessing()
      real(dp) :: x(12, 1), y(12)
      integer :: ncat(1), i, status
      type(ranger_options) :: options
      type(ranger_regression_forest) :: forest

      ncat = 3
      do i = 1, 4
         x(i, 1) = 1.0_dp
         y(i) = 10.0_dp + 0.1_dp * real(i, dp)
      end do
      do i = 5, 8
         x(i, 1) = 2.0_dp
         y(i) = 0.1_dp * real(i - 4, dp)
      end do
      do i = 9, 12
         x(i, 1) = 3.0_dp
         y(i) = 5.0_dp + 0.1_dp * real(i - 8, dp)
      end do
      options%num_trees = 1
      options%mtry = 1
      options%replace = .false.
      options%sample_fraction = 1.0_dp
      options%min_node_size = 1
      options%respect_unordered_factors = RANGER_UNORDERED_ORDER
      options%seed = 1_i64
      call fit_ranger_regression(x, y, forest, ncat=ncat, options=options, status=status)
      call require(status == 0, 'ordered-factor regression fit')
      call require(all(forest%category_map(:, 1) == [3, 1, 2]), 'global response-mean factor order')
   end subroutine test_ordered_factor_preprocessing

   subroutine test_survival_factor_ordering()
      real(dp) :: x(4, 1), time(4)
      real(dp), allocatable :: transformed(:,:)
      integer :: event(4), ncat(1)
      integer, allocatable :: mapping(:,:)

      ! Category 1 has survival exactly 0.5 from time 1 through the final
      ! censoring time 100, so survival::quantile.survfit returns 50.5.
      ! Category 2 is exactly 0.5 from 10 to 20, so it returns 15.  The
      ! midpoint rule therefore orders category 2 before category 1.
      x(:, 1) = [1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp]
      time = [1.0_dp, 100.0_dp, 10.0_dp, 20.0_dp]
      event = [1, 0, 1, 1]
      ncat = 2
      call order_factors_survival(x, time, event, ncat, mapping, transformed)
      call require(all(mapping(:, 1) == [2, 1]), 'survival factor flat-quantile midpoint order')
      call require(all(nint(transformed(:, 1)) == [2, 2, 1, 1]), 'survival factor mapping applied')
   end subroutine test_survival_factor_ordering

   subroutine test_survival_time_interest_parity()
      real(dp) :: x(8, 1), time(8), explicit_grid(4)
      integer :: event(8), i, status
      type(ranger_options) :: options
      type(ranger_survival_forest) :: forest

      do i = 1, 8
         x(i, 1) = real(i, dp)
      end do
      time = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      event = 1
      options%num_trees = 1
      options%mtry = 1
      options%replace = .false.
      options%sample_fraction = 0.5_dp
      options%min_node_size = 1
      options%min_bucket = 1
      options%seed = 4321_i64

      explicit_grid = [3.0_dp, 1.0_dp, 2.0_dp, 2.0_dp]
      call fit_ranger_survival(x, time, event, forest, options=options, time_interest=explicit_grid, status=status)
      call require(status == 0, 'survival time.interest explicit grid fit')
      call require(forest%ntime == 3, 'survival time.interest explicit deduplication')
      call require(maxval(abs(forest%unique_timepoints - [1.0_dp, 2.0_dp, 3.0_dp])) < 1.0e-15_dp, &
         'survival time.interest explicit sort/unique')

      call fit_ranger_survival(x, time, event, forest, options=options, time_interest=3, status=status)
      call require(status == 0, 'survival time.interest integer count fit')
      call require(forest%ntime == 3, 'survival time.interest integer count')
      call require(maxval(abs(forest%unique_timepoints - [1.0_dp, 2.0_dp, 4.0_dp])) < 1.0e-15_dp, &
         'survival time.interest R round-to-even count grid')

      call fit_ranger_survival(x, time, event, forest, options=options, time_interest=3.0_dp, status=status)
      call require(status == 0, 'survival time.interest real count fit')
      call require(maxval(abs(forest%unique_timepoints - [1.0_dp, 2.0_dp, 4.0_dp])) < 1.0e-15_dp, &
         'survival time.interest real count grid')
   end subroutine test_survival_time_interest_parity

   subroutine test_oob_quantiles()
      real(dp) :: x(80, 2), y(80), q(80, 3), probabilities(3)
      integer :: i, status
      type(ranger_options) :: options
      type(ranger_regression_forest) :: forest

      do i = 1, 80
         x(i, 1) = real(i, dp) / 80.0_dp
         x(i, 2) = sin(real(i, dp))
         y(i) = 2.0_dp * x(i, 1) + 0.1_dp * x(i, 2)
      end do
      options%num_trees = 200
      options%mtry = 2
      options%min_node_size = 5
      options%min_bucket = 5
      options%quantreg = .true.
      options%keep_inbag = .true.
      options%seed = 54321_i64
      call fit_ranger_regression(x, y, forest, options=options, status=status)
      call require(status == 0, 'OOB quantile preparation')
      call require(forest%quantile_oob_count >= 10, 'OOB quantile minimum tree count')
      probabilities = [0.1_dp, 0.5_dp, 0.9_dp]
      call predict_ranger_quantiles_oob(forest, probabilities, q)
      call require(all(q(:, 1) <= q(:, 2)) .and. all(q(:, 2) <= q(:, 3)), 'OOB quantile ordering')
   end subroutine test_oob_quantiles

   subroutine test_survival()
      real(dp) :: x(30, 1), time(30), chf(30, 10), surv(30, 10), risk_early, risk_late
      real(dp) :: grid(10)
      integer :: event(30), i, status
      type(ranger_options) :: options
      type(ranger_survival_forest) :: forest

      do i = 1, 15
         x(i, 1) = 0.0_dp
         time(i) = 1.0_dp + 0.2_dp * real(i - 1, dp)
         event(i) = 1
      end do
      do i = 16, 30
         x(i, 1) = 1.0_dp
         time(i) = 6.0_dp + 0.2_dp * real(i - 16, dp)
         event(i) = merge(1, 0, mod(i, 3) /= 0)
      end do
      grid = [(real(i, dp), i = 1, 10)]
      options%num_trees = 25
      options%mtry = 1
      options%min_node_size = 3
      options%min_bucket = 2
      options%seed = 2468
      call fit_ranger_survival(x, time, event, forest, options=options, time_interest=grid, status=status)
      call require(status == 0, 'survival fit status')
      call predict_ranger_survival(forest, x, chf, surv)
      risk_early = sum(chf(1:15, :)) / 15.0_dp
      risk_late = sum(chf(16:30, :)) / 15.0_dp
      call require(risk_early > risk_late, 'survival risk ordering')
      call require(all(surv >= 0.0_dp) .and. all(surv <= 1.0_dp), 'survival probabilities')
   end subroutine test_survival

   subroutine test_survival_concordance_parity()
      real(dp) :: time(3), risk(3), case_error(3), cindex
      integer :: status(3)

      ! The first two observations have the same time but different status.
      ! ranger treats that pair as permissible; with tied risk it contributes
      ! 0.5.  The other two permissible pairs are concordant.
      time = [1.0_dp, 1.0_dp, 2.0_dp]
      status = [1, 0, 1]
      risk = [3.0_dp, 3.0_dp, 1.0_dp]
      cindex = concordance_index(time, status, risk)
      call require(abs(cindex - 0.75_dp) < 1.0e-15_dp, 'survival concordance tied-time semantics')
      call concordance_casewise(time, status, risk, case_error)
      call require(abs(case_error(1) - 0.25_dp) < 1.0e-15_dp, 'survival casewise concordance event')
      call require(abs(case_error(2) - 0.5_dp) < 1.0e-15_dp, 'survival casewise concordance censor')
      call require(abs(case_error(3)) < 1.0e-15_dp, 'survival casewise concordance later event')
   end subroutine test_survival_concordance_parity

   subroutine test_utilities()
      real(dp) :: importance(5), pvalue(5), null_imp(5, 3)
      real(dp) :: pred(2, 4), mean_pred(2), var_pred(2), weights(3)
      integer :: inbag(3, 4), train_terminal(3, 2), test_terminal(2), status
      integer :: sizes(2), usage(2)
      type(ranger_options) :: options
      type(ranger_regression_forest) :: forest
      real(dp) :: x(12, 2), y(12)
      integer :: i

      importance = [-0.3_dp, -0.1_dp, 0.0_dp, 0.2_dp, 0.8_dp]
      call importance_pvalues_janitza(importance, pvalue, status)
      call require(status == 0, 'Janitza p-value status')
      call require(maxval(abs(pvalue - [1.0_dp, 0.8_dp, 0.6_dp, 0.2_dp, 0.0_dp])) < 1.0e-14_dp, &
         'Janitza strict-lower-bound p-values')
      null_imp(:, 1) = [-0.2_dp, 0.0_dp, 0.1_dp, 0.1_dp, 0.3_dp]
      null_imp(:, 2) = [-0.1_dp, 0.1_dp, 0.0_dp, 0.3_dp, 0.5_dp]
      null_imp(:, 3) = [0.0_dp, -0.1_dp, 0.2_dp, 0.2_dp, 0.4_dp]
      call importance_pvalues_altmann(importance, null_imp, pvalue)
      call require(all(pvalue > 0.0_dp) .and. all(pvalue <= 1.0_dp), 'Altmann p-value range')

      pred = reshape([1.0_dp, 2.0_dp, 1.2_dp, 2.1_dp, 0.9_dp, 1.8_dp, 1.1_dp, 2.2_dp], [2, 4])
      inbag = reshape([1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1], [3, 4])
      call infinitesimal_jackknife(pred, inbag, mean_pred, var_pred, calibrate=.true.)
      call require(all(ieee_is_finite(mean_pred)) .and. all(var_pred >= 0.0_dp), 'infinitesimal jackknife')

      train_terminal = reshape([1, 2, 1, 3, 2, 3], [3, 2])
      test_terminal = [1, 3]
      call case_specific_weights(train_terminal, test_terminal, weights)
      call require(abs(sum(weights) - 1.0_dp) < 1.0e-12_dp, 'case-specific weights normalization')

      do i = 1, 12
         x(i, 1) = real(i, dp)
         x(i, 2) = real(mod(i, 2), dp)
         y(i) = 0.5_dp * x(i, 1)
      end do
      options%num_trees = 2
      options%mtry = 2
      options%replace = .false.
      options%sample_fraction = 1.0_dp
      options%min_node_size = 2
      call fit_ranger_regression(x, y, forest, options=options, status=status)
      call require(status == 0, 'utility forest fit')
      call tree_sizes(forest%trees, sizes)
      call variable_usage(forest%trees, 2, usage)
      call require(all(sizes > 0), 'tree sizes')
      call require(sum(usage) > 0, 'variable usage')
   end subroutine test_utilities


   subroutine test_ij_small_sample()
      integer, parameter :: ntest = 20, ntrain = 6, ntree = 12
      real(dp) :: pred(ntest, ntree), raw(ntest), requested(ntest), mean1(ntest), mean2(ntest)
      integer :: inbag(ntrain, ntree), i, j, k

      do j = 1, ntree
         do i = 1, ntest
            pred(i, j) = cos(0.13_dp * real(i * j, dp)) + 0.01_dp * real(i, dp)
         end do
         inbag(:, j) = 0
         do k = 1, ntrain
            if (mod(j + 2 * k, 3) == 0) inbag(k, j) = 2
         end do
      end do
      call infinitesimal_jackknife(pred, inbag, mean1, raw, calibrate=.false.)
      call infinitesimal_jackknife(pred, inbag, mean2, requested, calibrate=.true., seed=19_i64)
      call require(maxval(abs(mean1 - mean2)) < 1.0e-14_dp, 'IJ small-sample mean unchanged')
      call require(maxval(abs(raw - requested)) < 1.0e-14_dp, 'IJ calibration disabled for <=20 predictions')
   end subroutine test_ij_small_sample

   subroutine test_ij_calibration()
      integer, parameter :: ntest = 21, ntrain = 8, ntree = 48
      real(dp) :: pred(ntest, ntree), raw(ntest), calibrated(ntest), mean1(ntest), mean2(ntest)
      integer :: inbag(ntrain, ntree), i, j, k

      do j = 1, ntree
         do i = 1, ntest
            pred(i, j) = 0.03_dp * real(i, dp) + sin(0.17_dp * real(i * j, dp))
         end do
         inbag(:, j) = 0
         do k = 1, ntrain
            inbag(k, j) = merge(2, 0, mod(j + 3 * k, ntrain) < 4)
         end do
      end do
      call infinitesimal_jackknife(pred, inbag, mean1, raw, calibrate=.false.)
      call infinitesimal_jackknife(pred, inbag, mean2, calibrated, calibrate=.true., seed=12345_i64)
      call require(maxval(abs(mean1 - mean2)) < 1.0e-14_dp, 'IJ calibrated mean unchanged')
      call require(all(ieee_is_finite(calibrated)) .and. all(calibrated >= 0.0_dp), 'IJ EB calibration finite')
      call require(maxval(abs(calibrated - raw)) > 1.0e-4_dp, 'IJ EB calibration changes raw variance')
   end subroutine test_ij_calibration

   integer function first_terminal(status, nnode) result(node)
      integer, intent(in) :: status(:), nnode
      integer :: i
      node = 1
      do i = 1, nnode
         if (status(i) /= 1) then
            node = i
            return
         end if
      end do
   end function first_terminal

   subroutine require(condition, name)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      if (.not. condition) then
         print '(a)', 'FAILED: ' // trim(name)
         error stop 1
      end if
   end subroutine require

end program test_ranger
