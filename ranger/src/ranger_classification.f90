! SPDX-License-Identifier: GPL-3.0-only
! Modern Fortran translation of ranger 0.18.0 computational code.
! Upstream C++ core: Copyright (c) 2014-2018 Marvin N. Wright; MIT licensed.
! The upstream R package is GPL-3. See NOTICE.md for full provenance.
module ranger_classification
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp, i64
   use ranger_rng, only : ranger_rng_state, shuffle_int, shuffle_real
   use ranger_types, only : ranger_options, ranger_classification_forest, ranger_probability_forest, ranger_tree
   use ranger_types, only : RANGER_IMPORTANCE_IMPURITY, RANGER_IMPORTANCE_IMPURITY_CORRECTED
   use ranger_types, only : RANGER_IMPORTANCE_PERMUTATION, RANGER_IMPORTANCE_PERMUTATION_CASEWISE
   use ranger_types, only : RANGER_SPLIT_STANDARD, RANGER_SPLIT_EXTRATREES, RANGER_SPLIT_HELLINGER
   use ranger_types, only : RANGER_UNORDERED_IGNORE, RANGER_UNORDERED_ORDER, RANGER_UNORDERED_PARTITION
   use ranger_factor_order, only : order_factors_classification, identity_factor_map, apply_factor_map
   use ranger_forest_common, only : prepare_categories, make_case_weights, draw_inbag, draw_classwise_inbag
   use ranger_forest_common, only : inbag_to_observations, argmax_real, set_ok, set_error, validate_options
   use ranger_forest_common, only : prepare_training_predictors
   use ranger_forest_common, only : initialize_forest_rng, initialize_tree_rng, tree_options_for_index
   use ranger_tree_classification, only : build_classification_tree, predict_classification_tree
   use ranger_tree_common, only : tree_uses_variable
   implicit none
   private

   public :: fit_ranger_classification, predict_ranger_classification
   public :: fit_ranger_probability, predict_ranger_probability

contains

   subroutine fit_ranger_classification(x, y, forest, ncat, options, class_weights, case_weights, &
      sample_fraction_class, min_node_size_class, min_bucket_class, preset_inbag, status, message)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: y(:)
      type(ranger_classification_forest), intent(out) :: forest
      integer, intent(in), optional :: ncat(:)
      type(ranger_options), intent(in), optional :: options
      real(dp), intent(in), optional :: class_weights(:), case_weights(:), sample_fraction_class(:)
      integer, intent(in), optional :: min_node_size_class(:), min_bucket_class(:), preset_inbag(:,:)
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      type(ranger_options) :: opt, tree_opt
      type(ranger_rng_state) :: rng, tree_rng
      integer, allocatable :: cats(:), cats_train(:), category_map(:,:), inbag(:), obs_index(:), tree_pred(:), leaf(:)
      integer, allocatable :: min_node(:), min_bucket(:)
      real(dp), allocatable :: obs_weight(:), cweight(:), case_weight(:), imp_sum(:), imp_sumsq(:), x_train(:,:), x_ordered(:,:)
      real(dp), allocatable :: local_sum(:,:)
      logical, allocatable :: used_global(:)
      integer :: n, p, nclass, t, i, k, rc, noob, nerr
      logical :: use_case_weights
      character(len=200) :: message_local

      call set_ok(status, message)
      n = size(x, 1)
      p = size(x, 2)
      if (n <= 1 .or. p <= 0 .or. size(y) /= n) then
         call set_error(1, 'x and y have incompatible or empty dimensions', status, message)
         return
      end if
      if (any((.not. ieee_is_finite(x)) .and. (.not. ieee_is_nan(x)))) then
         call set_error(2, 'x contains infinite values', status, message)
         return
      end if
      if (minval(y) /= 1) then
         call set_error(3, 'class labels must be consecutive positive integers beginning at 1', status, message)
         return
      end if
      nclass = maxval(y)
      if (nclass < 2) then
         call set_error(4, 'classification requires at least two classes', status, message)
         return
      end if
      do k = 1, nclass
         if (count(y == k) == 0) then
            call set_error(5, 'class labels must contain every integer through max(y)', status, message)
            return
         end if
      end do

      opt = ranger_options()
      if (present(options)) opt = options
      if (opt%mtry <= 0) opt%mtry = max(1, int(sqrt(real(p, dp))))
      call validate_options(opt, p, rc, message_local)
      if (rc /= 0) then
         call set_error(10 + rc, trim(message_local), status, message)
         return
      end if
      if (opt%split_rule /= RANGER_SPLIT_STANDARD .and. opt%split_rule /= RANGER_SPLIT_EXTRATREES .and. &
         opt%split_rule /= RANGER_SPLIT_HELLINGER) then
         call set_error(18, 'classification supports standard, extratrees, or hellinger split rules', status, message)
         return
      end if
      if (opt%split_rule == RANGER_SPLIT_HELLINGER .and. nclass /= 2) then
         call set_error(19, 'hellinger splitrule requires exactly two classes', status, message)
         return
      end if
      call prepare_categories(x, ncat, opt%na_learn, cats, rc)
      if (rc /= 0) then
         call set_error(20 + rc, 'invalid categorical predictor coding or missing-value policy', status, message)
         return
      end if
      if (opt%respect_unordered_factors == RANGER_UNORDERED_PARTITION .and. &
         any(cats > digits(1.0_dp))) then
         call set_error(28, 'unordered partition factors are limited to 53 observed levels', status, message)
         return
      end if
      if (opt%respect_unordered_factors == RANGER_UNORDERED_ORDER) then
         call order_factors_classification(x, y, nclass, cats, category_map, x_ordered)
         opt%respect_unordered_factors = RANGER_UNORDERED_IGNORE
      else
         call identity_factor_map(cats, category_map)
         call apply_factor_map(x, cats, category_map, x_ordered)
      end if
      call make_case_weights(n, case_weights, case_weight, rc, use_case_weights)
      if (rc /= 0) then
         call set_error(30 + rc, 'case_weights must be nonnegative with positive total weight', status, message)
         return
      end if
      if (opt%holdout .and. .not. use_case_weights) then
         call set_error(33, 'nonconstant case_weights are required in holdout mode', status, message)
         return
      end if
      call prepare_class_controls(nclass, opt, class_weights, min_node_size_class, min_bucket_class, cweight, &
         min_node, min_bucket, rc)
      if (rc /= 0) then
         call set_error(40 + rc, 'invalid class weights or class-wise node-size controls', status, message)
         return
      end if
      if (present(sample_fraction_class)) then
         if (size(sample_fraction_class) /= nclass) then
            call set_error(50, 'sample_fraction_class must have one entry per class', status, message)
            return
         end if
         if (present(case_weights)) then
            call set_error(52, 'case_weights and class-wise sampling cannot be combined', status, message)
            return
         end if
         if (.not. opt%replace) then
            do k = 1, nclass
               if (sample_fraction_class(k) * real(n, dp) > real(count(y == k), dp)) then
                  call set_error(53, 'class-wise sample fraction requests more observations than available', status, message)
                  return
               end if
            end do
         end if
      end if
      if (present(preset_inbag)) then
         if (size(preset_inbag, 1) /= n .or. size(preset_inbag, 2) /= opt%num_trees .or. any(preset_inbag < 0)) then
            call set_error(51, 'preset_inbag must have shape (nobs,num_trees) and be nonnegative', status, message)
            return
         end if
      end if
      if (present(preset_inbag) .and. use_case_weights) then
         call set_error(54, 'case_weights and preset_inbag cannot be combined', status, message)
         return
      end if
      if (present(preset_inbag) .and. present(sample_fraction_class)) then
         call set_error(55, 'class-wise sampling and preset_inbag cannot be combined', status, message)
         return
      end if

      call initialize_classification_forest(forest, n, p, nclass, cats, opt)
      forest%category_map = category_map
      allocate(inbag(n), tree_pred(n), leaf(n), used_global(p), imp_sum(p), imp_sumsq(p))
      imp_sum = 0.0_dp
      imp_sumsq = 0.0_dp
      used_global = .false.
      if (opt%local_importance .or. opt%importance_mode == RANGER_IMPORTANCE_PERMUTATION_CASEWISE) then
         allocate(local_sum(p, n))
         local_sum = 0.0_dp
      end if
      call initialize_forest_rng(rng, opt%seed)
      call prepare_training_predictors(rng, x_ordered, cats, opt, x_train, cats_train)

      do t = 1, opt%num_trees
         call initialize_tree_rng(rng, tree_rng, opt%seed, t)
         if (use_case_weights) then
            call draw_inbag(tree_rng, opt, case_weight, .true., inbag, rc)
         else if (present(sample_fraction_class)) then
            call draw_classwise_inbag(tree_rng, y, nclass, opt, sample_fraction_class, inbag, rc)
         else if (present(preset_inbag)) then
            inbag = preset_inbag(:, t)
            rc = 0
         else
            call draw_inbag(tree_rng, opt, case_weight, .false., inbag, rc)
         end if
         if (rc /= 0) then
            call set_error(60, 'sampling failed', status, message)
            return
         end if
         if (allocated(forest%inbag)) forest%inbag(:, t) = inbag
         call inbag_to_observations(inbag, obs_index, obs_weight)
         if (present(preset_inbag) .and. .not. use_case_weights .and. .not. present(sample_fraction_class)) &
            call shuffle_int(tree_rng, obs_index)
         call tree_options_for_index(opt, t, tree_opt)
         call build_classification_tree(x_train, y, obs_index, obs_weight, cats_train, nclass, cweight, tree_opt, &
            min_node, min_bucket, &
            tree_rng, used_global, forest%trees(t))
         if (opt%oob_error) then
            call predict_classification_tree(forest%trees(t), x_ordered, cats, tree_pred, terminal_node=leaf)
            do i = 1, n
               if (.not. is_evaluation_oob(i, inbag, case_weight, opt%holdout)) cycle
               forest%oob_votes(i, tree_pred(i)) = forest%oob_votes(i, tree_pred(i)) + 1.0_dp
               forest%oob_count(i) = forest%oob_count(i) + 1
            end do
            call update_classification_error(y, forest%oob_votes, forest%oob_count, forest%oob_prediction, &
               forest%prediction_error(t), rng)
         end if
         if (opt%importance_mode == RANGER_IMPORTANCE_IMPURITY .or. &
            opt%importance_mode == RANGER_IMPORTANCE_IMPURITY_CORRECTED) then
            imp_sum = imp_sum + forest%trees(t)%impurity_decrease
         else if (opt%importance_mode == RANGER_IMPORTANCE_PERMUTATION .or. &
            opt%importance_mode == RANGER_IMPORTANCE_PERMUTATION_CASEWISE) then
            if (allocated(local_sum)) then
               call classification_permutation_importance(forest%trees(t), x_ordered, y, cats, inbag, case_weight, opt, tree_rng, &
                  imp_sum, imp_sumsq, local_sum)
            else
               call classification_permutation_importance(forest%trees(t), x_ordered, y, cats, inbag, case_weight, opt, tree_rng, &
                  imp_sum, imp_sumsq)
            end if
         end if
         deallocate(obs_index, obs_weight)
      end do

      if (opt%importance_mode == RANGER_IMPORTANCE_IMPURITY .or. &
         opt%importance_mode == RANGER_IMPORTANCE_IMPURITY_CORRECTED) then
         forest%variable_importance = imp_sum / real(opt%num_trees, dp)
      else if (opt%importance_mode == RANGER_IMPORTANCE_PERMUTATION .or. &
         opt%importance_mode == RANGER_IMPORTANCE_PERMUTATION_CASEWISE) then
         call finish_importance(imp_sum, imp_sumsq, opt%num_trees, opt%scale_permutation_importance, &
            forest%variable_importance, forest%variable_importance_sd)
      end if
      if (allocated(local_sum)) forest%local_importance = local_sum / real(opt%num_trees, dp)

      if (opt%oob_error) then
         noob = count(forest%oob_count > 0)
         nerr = count(forest%oob_count > 0 .and. forest%oob_prediction /= y)
         if (noob > 0 .and. size(forest%prediction_error) > 0) then
            forest%prediction_error(size(forest%prediction_error)) = real(nerr, dp) / real(noob, dp)
         end if
      end if

   end subroutine fit_ranger_classification

   subroutine predict_ranger_classification(forest, x, prediction, vote_fraction, per_tree, terminal_nodes, seed)
      type(ranger_classification_forest), intent(in) :: forest
      real(dp), intent(in) :: x(:,:)
      integer, intent(out) :: prediction(:)
      real(dp), intent(out), optional :: vote_fraction(:,:)
      integer, intent(out), optional :: per_tree(:,:), terminal_nodes(:,:)
      integer(i64), intent(in), optional :: seed
      integer, allocatable :: tree_pred(:), leaf(:)
      real(dp), allocatable :: votes(:,:), x_mapped(:,:)
      type(ranger_rng_state) :: tie_rng
      integer(i64) :: tie_seed
      integer :: n, t, i

      n = size(x, 1)
      if (size(x, 2) /= forest%nvar .or. size(prediction) /= n) &
         error stop 'predict_ranger_classification: incompatible dimensions'
      if (present(vote_fraction)) then
         if (size(vote_fraction, 1) /= n .or. size(vote_fraction, 2) /= forest%nclass) &
            error stop 'predict_ranger_classification: vote_fraction has wrong shape'
      end if
      if (present(per_tree)) then
         if (size(per_tree, 1) /= n .or. size(per_tree, 2) /= size(forest%trees)) &
            error stop 'predict_ranger_classification: per_tree has wrong shape'
      end if
      if (present(terminal_nodes)) then
         if (size(terminal_nodes, 1) /= n .or. size(terminal_nodes, 2) /= size(forest%trees)) &
            error stop 'predict_ranger_classification: terminal_nodes has wrong shape'
      end if
      allocate(tree_pred(n), leaf(n), votes(n, forest%nclass))
      call apply_factor_map(x, forest%ncat, forest%category_map, x_mapped)
      tie_seed = 0_i64
      if (present(seed)) tie_seed = seed
      call initialize_forest_rng(tie_rng, tie_seed)
      votes = 0.0_dp
      do t = 1, size(forest%trees)
         call predict_classification_tree(forest%trees(t), x_mapped, forest%ncat, tree_pred, terminal_node=leaf)
         do i = 1, n
            votes(i, tree_pred(i)) = votes(i, tree_pred(i)) + 1.0_dp
         end do
         if (present(per_tree)) per_tree(:, t) = tree_pred
         if (present(terminal_nodes)) terminal_nodes(:, t) = leaf
      end do
      do i = 1, n
         prediction(i) = random_argmax(votes(i, :), tie_rng)
      end do
      if (present(vote_fraction)) vote_fraction = votes / real(size(forest%trees), dp)
   end subroutine predict_ranger_classification

   subroutine fit_ranger_probability(x, y, forest, ncat, options, class_weights, case_weights, &
      sample_fraction_class, min_node_size_class, min_bucket_class, preset_inbag, status, message)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: y(:)
      type(ranger_probability_forest), intent(out) :: forest
      integer, intent(in), optional :: ncat(:)
      type(ranger_options), intent(in), optional :: options
      real(dp), intent(in), optional :: class_weights(:), case_weights(:), sample_fraction_class(:)
      integer, intent(in), optional :: min_node_size_class(:), min_bucket_class(:), preset_inbag(:,:)
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      type(ranger_options) :: opt, tree_opt
      type(ranger_rng_state) :: rng, tree_rng
      integer, allocatable :: cats(:), cats_train(:), category_map(:,:), inbag(:), obs_index(:), leaf(:), min_node(:), min_bucket(:)
      real(dp), allocatable :: obs_weight(:), cweight(:), case_weight(:), tree_prob(:,:), imp_sum(:), imp_sumsq(:)
      real(dp), allocatable :: x_train(:,:), x_ordered(:,:)
      real(dp), allocatable :: local_sum(:,:)
      logical, allocatable :: used_global(:)
      integer :: n, p, nclass, t, i, k, rc
      logical :: use_case_weights
      character(len=200) :: message_local

      call set_ok(status, message)
      n = size(x, 1)
      p = size(x, 2)
      if (n <= 1 .or. p <= 0 .or. size(y) /= n) then
         call set_error(1, 'x and y have incompatible or empty dimensions', status, message)
         return
      end if
      if (any((.not. ieee_is_finite(x)) .and. (.not. ieee_is_nan(x)))) then
         call set_error(2, 'x contains infinite values', status, message)
         return
      end if
      if (minval(y) /= 1) then
         call set_error(3, 'class labels must be consecutive positive integers beginning at 1', status, message)
         return
      end if
      nclass = maxval(y)
      if (nclass < 2) then
         call set_error(4, 'probability forests require at least two classes', status, message)
         return
      end if
      do k = 1, nclass
         if (count(y == k) == 0) then
            call set_error(5, 'class labels must contain every integer through max(y)', status, message)
            return
         end if
      end do
      opt = ranger_options()
      if (present(options)) opt = options
      if (opt%mtry <= 0) opt%mtry = max(1, int(sqrt(real(p, dp))))
      call validate_options(opt, p, rc, message_local)
      if (rc /= 0) then
         call set_error(10 + rc, trim(message_local), status, message)
         return
      end if
      if (opt%split_rule /= RANGER_SPLIT_STANDARD .and. opt%split_rule /= RANGER_SPLIT_EXTRATREES .and. &
         opt%split_rule /= RANGER_SPLIT_HELLINGER) then
         call set_error(18, 'classification supports standard, extratrees, or hellinger split rules', status, message)
         return
      end if
      if (opt%split_rule == RANGER_SPLIT_HELLINGER .and. nclass /= 2) then
         call set_error(19, 'hellinger splitrule requires exactly two classes', status, message)
         return
      end if
      call prepare_categories(x, ncat, opt%na_learn, cats, rc)
      if (rc /= 0) then
         call set_error(20 + rc, 'invalid categorical predictor coding or missing-value policy', status, message)
         return
      end if
      if (opt%respect_unordered_factors == RANGER_UNORDERED_PARTITION .and. &
         any(cats > digits(1.0_dp))) then
         call set_error(28, 'unordered partition factors are limited to 53 observed levels', status, message)
         return
      end if
      if (opt%respect_unordered_factors == RANGER_UNORDERED_ORDER) then
         call order_factors_classification(x, y, nclass, cats, category_map, x_ordered)
         opt%respect_unordered_factors = RANGER_UNORDERED_IGNORE
      else
         call identity_factor_map(cats, category_map)
         call apply_factor_map(x, cats, category_map, x_ordered)
      end if
      call make_case_weights(n, case_weights, case_weight, rc, use_case_weights)
      if (rc /= 0) then
         call set_error(30 + rc, 'case_weights must be nonnegative with positive total weight', status, message)
         return
      end if
      if (opt%holdout .and. .not. use_case_weights) then
         call set_error(33, 'nonconstant case_weights are required in holdout mode', status, message)
         return
      end if
      call prepare_class_controls(nclass, opt, class_weights, min_node_size_class, min_bucket_class, cweight, &
         min_node, min_bucket, rc, probability_default=.true.)
      if (rc /= 0) then
         call set_error(40 + rc, 'invalid class weights or class-wise node-size controls', status, message)
         return
      end if
      if (present(sample_fraction_class)) then
         if (size(sample_fraction_class) /= nclass) then
            call set_error(50, 'sample_fraction_class must have one entry per class', status, message)
            return
         end if
         if (present(case_weights)) then
            call set_error(52, 'case_weights and class-wise sampling cannot be combined', status, message)
            return
         end if
         if (.not. opt%replace) then
            do k = 1, nclass
               if (sample_fraction_class(k) * real(n, dp) > real(count(y == k), dp)) then
                  call set_error(53, 'class-wise sample fraction requests more observations than available', status, message)
                  return
               end if
            end do
         end if
      end if
      if (present(preset_inbag)) then
         if (size(preset_inbag, 1) /= n .or. size(preset_inbag, 2) /= opt%num_trees .or. any(preset_inbag < 0)) then
            call set_error(51, 'preset_inbag must have shape (nobs,num_trees) and be nonnegative', status, message)
            return
         end if
      end if
      if (present(preset_inbag) .and. use_case_weights) then
         call set_error(54, 'case_weights and preset_inbag cannot be combined', status, message)
         return
      end if
      if (present(preset_inbag) .and. present(sample_fraction_class)) then
         call set_error(55, 'class-wise sampling and preset_inbag cannot be combined', status, message)
         return
      end if

      call initialize_probability_forest(forest, n, p, nclass, cats, opt)
      forest%category_map = category_map
      allocate(inbag(n), leaf(n), tree_prob(n, nclass), used_global(p), imp_sum(p), imp_sumsq(p))
      imp_sum = 0.0_dp
      imp_sumsq = 0.0_dp
      used_global = .false.
      if (opt%local_importance .or. opt%importance_mode == RANGER_IMPORTANCE_PERMUTATION_CASEWISE) then
         allocate(local_sum(p, n))
         local_sum = 0.0_dp
      end if
      call initialize_forest_rng(rng, opt%seed)
      call prepare_training_predictors(rng, x_ordered, cats, opt, x_train, cats_train)
      do t = 1, opt%num_trees
         call initialize_tree_rng(rng, tree_rng, opt%seed, t)
         if (use_case_weights) then
            call draw_inbag(tree_rng, opt, case_weight, .true., inbag, rc)
         else if (present(sample_fraction_class)) then
            call draw_classwise_inbag(tree_rng, y, nclass, opt, sample_fraction_class, inbag, rc)
         else if (present(preset_inbag)) then
            inbag = preset_inbag(:, t)
            rc = 0
         else
            call draw_inbag(tree_rng, opt, case_weight, .false., inbag, rc)
         end if
         if (rc /= 0) then
            call set_error(60, 'sampling failed', status, message)
            return
         end if
         if (allocated(forest%inbag)) forest%inbag(:, t) = inbag
         call inbag_to_observations(inbag, obs_index, obs_weight)
         if (present(preset_inbag) .and. .not. use_case_weights .and. .not. present(sample_fraction_class)) &
            call shuffle_int(tree_rng, obs_index)
         call tree_options_for_index(opt, t, tree_opt)
         call build_classification_tree(x_train, y, obs_index, obs_weight, cats_train, nclass, cweight, tree_opt, &
            min_node, min_bucket, &
            tree_rng, used_global, forest%trees(t))
         if (opt%oob_error) then
            call predict_classification_tree(forest%trees(t), x_ordered, cats, probability=tree_prob, terminal_node=leaf)
            do i = 1, n
               if (.not. is_evaluation_oob(i, inbag, case_weight, opt%holdout)) cycle
               forest%oob_probability(i, :) = forest%oob_probability(i, :) + tree_prob(i, :)
               forest%oob_count(i) = forest%oob_count(i) + 1
            end do
            call update_probability_error(y, forest%oob_probability, forest%oob_count, forest%prediction_error(t))
         end if
         if (opt%importance_mode == RANGER_IMPORTANCE_IMPURITY .or. &
            opt%importance_mode == RANGER_IMPORTANCE_IMPURITY_CORRECTED) then
            imp_sum = imp_sum + forest%trees(t)%impurity_decrease
         else if (opt%importance_mode == RANGER_IMPORTANCE_PERMUTATION .or. &
            opt%importance_mode == RANGER_IMPORTANCE_PERMUTATION_CASEWISE) then
            if (allocated(local_sum)) then
               call probability_permutation_importance(forest%trees(t), x_ordered, y, cats, inbag, case_weight, opt, tree_rng, &
                  imp_sum, imp_sumsq, local_sum)
            else
               call probability_permutation_importance(forest%trees(t), x_ordered, y, cats, inbag, case_weight, opt, tree_rng, &
                  imp_sum, imp_sumsq)
            end if
         end if
         deallocate(obs_index, obs_weight)
      end do
      if (opt%oob_error) then
         do i = 1, n
            if (forest%oob_count(i) > 0) then
               forest%oob_probability(i, :) = forest%oob_probability(i, :) / real(forest%oob_count(i), dp)
            else
               forest%oob_probability(i, :) = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
         end do
      end if
      if (opt%importance_mode == RANGER_IMPORTANCE_IMPURITY .or. &
         opt%importance_mode == RANGER_IMPORTANCE_IMPURITY_CORRECTED) then
         forest%variable_importance = imp_sum / real(opt%num_trees, dp)
      else if (opt%importance_mode == RANGER_IMPORTANCE_PERMUTATION .or. &
         opt%importance_mode == RANGER_IMPORTANCE_PERMUTATION_CASEWISE) then
         call finish_importance(imp_sum, imp_sumsq, opt%num_trees, opt%scale_permutation_importance, &
            forest%variable_importance, forest%variable_importance_sd)
      end if
      if (allocated(local_sum)) forest%local_importance = local_sum / real(opt%num_trees, dp)
   end subroutine fit_ranger_probability

   subroutine predict_ranger_probability(forest, x, probability, prediction, per_tree_probability, terminal_nodes)
      type(ranger_probability_forest), intent(in) :: forest
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: probability(:,:)
      integer, intent(out), optional :: prediction(:), terminal_nodes(:,:)
      real(dp), intent(out), optional :: per_tree_probability(:,:,:)
      real(dp), allocatable :: tree_prob(:,:), x_mapped(:,:)
      integer, allocatable :: leaf(:)
      integer :: n, t, i

      n = size(x, 1)
      if (size(x, 2) /= forest%nvar .or. size(probability, 1) /= n .or. size(probability, 2) /= forest%nclass) &
         error stop 'predict_ranger_probability: incompatible dimensions'
      if (present(prediction)) then
         if (size(prediction) /= n) error stop 'predict_ranger_probability: prediction has wrong size'
      end if
      if (present(terminal_nodes)) then
         if (size(terminal_nodes, 1) /= n .or. size(terminal_nodes, 2) /= size(forest%trees)) &
            error stop 'predict_ranger_probability: terminal_nodes has wrong shape'
      end if
      if (present(per_tree_probability)) then
         if (size(per_tree_probability, 1) /= n .or. size(per_tree_probability, 2) /= forest%nclass .or. &
            size(per_tree_probability, 3) /= size(forest%trees)) &
            error stop 'predict_ranger_probability: per_tree_probability has wrong shape'
      end if
      allocate(tree_prob(n, forest%nclass), leaf(n))
      call apply_factor_map(x, forest%ncat, forest%category_map, x_mapped)
      probability = 0.0_dp
      do t = 1, size(forest%trees)
         call predict_classification_tree(forest%trees(t), x_mapped, forest%ncat, probability=tree_prob, terminal_node=leaf)
         probability = probability + tree_prob
         if (present(per_tree_probability)) per_tree_probability(:, :, t) = tree_prob
         if (present(terminal_nodes)) terminal_nodes(:, t) = leaf
      end do
      probability = probability / real(size(forest%trees), dp)
      if (present(prediction)) then
         do i = 1, n
            prediction(i) = argmax_real(probability(i, :))
         end do
      end if
   end subroutine predict_ranger_probability

   subroutine prepare_class_controls(nclass, options, class_weights, min_node_input, min_bucket_input, weights, &
      min_node, min_bucket, status, probability_default)
      integer, intent(in) :: nclass
      type(ranger_options), intent(in) :: options
      real(dp), intent(in), optional :: class_weights(:)
      integer, intent(in), optional :: min_node_input(:), min_bucket_input(:)
      real(dp), allocatable, intent(out) :: weights(:)
      integer, allocatable, intent(out) :: min_node(:), min_bucket(:)
      integer, intent(out) :: status
      logical, intent(in), optional :: probability_default
      integer :: node_default, bucket_default

      status = 0
      allocate(weights(nclass))
      weights = 1.0_dp
      if (present(class_weights)) then
         if (size(class_weights) /= nclass .or. any(class_weights < 0.0_dp) .or. sum(class_weights) <= 0.0_dp) then
            status = 1
            allocate(min_node(0), min_bucket(0))
            return
         end if
         weights = class_weights
      end if
      node_default = 1
      if (present(probability_default)) then
         if (probability_default) node_default = 10
      end if
      if (options%min_node_size > 0) node_default = options%min_node_size
      bucket_default = 1
      if (options%min_bucket > 0) bucket_default = options%min_bucket

      if (present(min_node_input)) then
         ! ranger treats a scalar zero as the "use default" sentinel, but a
         ! class-specific vector may contain genuine zero limits.
         if (size(min_node_input) /= nclass .or. any(min_node_input < 0)) then
            status = 2
            allocate(min_node(0), min_bucket(0))
            return
         end if
         allocate(min_node(nclass))
         min_node = min_node_input
      else
         allocate(min_node(1))
         min_node(1) = node_default
      end if
      if (present(min_bucket_input)) then
         ! As with min.node.size, zero is valid inside a class-specific vector.
         if (size(min_bucket_input) /= nclass .or. any(min_bucket_input < 0)) then
            status = 3
            allocate(min_bucket(0))
            return
         end if
         allocate(min_bucket(nclass))
         min_bucket = min_bucket_input
      else
         allocate(min_bucket(1))
         min_bucket(1) = bucket_default
      end if
   end subroutine prepare_class_controls

   subroutine initialize_classification_forest(forest, n, p, nclass, ncat, options)
      type(ranger_classification_forest), intent(out) :: forest
      integer, intent(in) :: n, p, nclass, ncat(:)
      type(ranger_options), intent(in) :: options
      forest%nclass = nclass
      forest%nvar = p
      forest%nobs = n
      allocate(forest%ncat(p), forest%category_map(max(1, maxval(ncat)), p))
      allocate(forest%trees(options%num_trees), forest%oob_prediction(n), forest%oob_count(n))
      allocate(forest%oob_votes(n, nclass), forest%prediction_error(options%num_trees))
      allocate(forest%variable_importance(p), forest%variable_importance_sd(p))
      forest%ncat = ncat
      forest%category_map = 0
      forest%oob_prediction = 0
      forest%oob_count = 0
      forest%oob_votes = 0.0_dp
      if (options%oob_error) then
         forest%prediction_error = 0.0_dp
      else
         forest%prediction_error = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
      forest%variable_importance = 0.0_dp
      forest%variable_importance_sd = 0.0_dp
      if (options%keep_inbag) then
         allocate(forest%inbag(n, options%num_trees))
         forest%inbag = 0
      end if
      if (options%local_importance .or. options%importance_mode == RANGER_IMPORTANCE_PERMUTATION_CASEWISE) then
         allocate(forest%local_importance(p, n))
         forest%local_importance = 0.0_dp
      end if
   end subroutine initialize_classification_forest

   subroutine initialize_probability_forest(forest, n, p, nclass, ncat, options)
      type(ranger_probability_forest), intent(out) :: forest
      integer, intent(in) :: n, p, nclass, ncat(:)
      type(ranger_options), intent(in) :: options
      forest%nclass = nclass
      forest%nvar = p
      forest%nobs = n
      allocate(forest%ncat(p), forest%category_map(max(1, maxval(ncat)), p))
      allocate(forest%trees(options%num_trees), forest%oob_probability(n, nclass), forest%oob_count(n))
      allocate(forest%prediction_error(options%num_trees), forest%variable_importance(p), forest%variable_importance_sd(p))
      forest%ncat = ncat
      forest%category_map = 0
      forest%oob_probability = 0.0_dp
      forest%oob_count = 0
      if (options%oob_error) then
         forest%prediction_error = 0.0_dp
      else
         forest%prediction_error = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
      forest%variable_importance = 0.0_dp
      forest%variable_importance_sd = 0.0_dp
      if (options%keep_inbag) then
         allocate(forest%inbag(n, options%num_trees))
         forest%inbag = 0
      end if
      if (options%local_importance .or. options%importance_mode == RANGER_IMPORTANCE_PERMUTATION_CASEWISE) then
         allocate(forest%local_importance(p, n))
         forest%local_importance = 0.0_dp
      end if
   end subroutine initialize_probability_forest

   logical function is_evaluation_oob(i, inbag, case_weight, holdout) result(use_case)
      integer, intent(in) :: i, inbag(:)
      real(dp), intent(in) :: case_weight(:)
      logical, intent(in) :: holdout
      if (holdout .and. any(case_weight <= 0.0_dp)) then
         use_case = case_weight(i) <= 0.0_dp
      else
         use_case = inbag(i) == 0
      end if
   end function is_evaluation_oob

   subroutine update_classification_error(y, votes, oob_count, prediction, error, rng)
      integer, intent(in) :: y(:), oob_count(:)
      real(dp), intent(in) :: votes(:,:)
      integer, intent(inout) :: prediction(:)
      real(dp), intent(out) :: error
      type(ranger_rng_state), intent(in) :: rng
      integer :: i, nuse, nerr
      nuse = 0
      nerr = 0
      error = ieee_value(0.0_dp, ieee_quiet_nan)
      do i = 1, size(y)
         if (oob_count(i) <= 0) cycle
         prediction(i) = random_argmax(votes(i, :), rng)
         nuse = nuse + 1
         if (prediction(i) /= y(i)) nerr = nerr + 1
      end do
      if (nuse > 0) error = real(nerr, dp) / real(nuse, dp)
   end subroutine update_classification_error

   integer function random_argmax(values, rng) result(index_value)
      real(dp), intent(in) :: values(:)
      type(ranger_rng_state), intent(in) :: rng
      type(ranger_rng_state) :: local_rng
      integer, allocatable :: tied(:)
      real(dp) :: maximum
      integer :: i, ntied

      maximum = maxval(values)
      allocate(tied(size(values)))
      ntied = 0
      do i = 1, size(values)
         if (abs(values(i) - maximum) <= 0.0_dp) then
            ntied = ntied + 1
            tied(ntied) = i
         end if
      end do
      if (ntied == 1) then
         index_value = tied(1)
      else
         ! ranger's mostFrequentValue() receives mt19937_64 by value, so tie
         ! resolution does not advance the forest RNG.  Class-container
         ! iteration order is a libstdc++ implementation detail; the native
         ! API uses ascending class IDs while preserving the exact uniform
         ! random choice among tied classes.
         local_rng = rng
         index_value = tied(local_rng%randint(ntied))
      end if
   end function random_argmax

   subroutine update_probability_error(y, prob_sum, oob_count, error)
      integer, intent(in) :: y(:), oob_count(:)
      real(dp), intent(in) :: prob_sum(:,:)
      real(dp), intent(out) :: error
      integer :: i, nuse
      real(dp) :: ptrue, error_sum
      error = ieee_value(0.0_dp, ieee_quiet_nan)
      error_sum = 0.0_dp
      nuse = 0
      do i = 1, size(y)
         if (oob_count(i) <= 0) cycle
         ptrue = prob_sum(i, y(i)) / real(oob_count(i), dp)
         error_sum = error_sum + (1.0_dp - ptrue) ** 2
         nuse = nuse + 1
      end do
      if (nuse > 0) error = error_sum / real(nuse, dp)
   end subroutine update_probability_error

   subroutine classification_permutation_importance(tree, x, y, ncat, inbag, case_weight, options, rng, &
      sum_imp, sumsq_imp, local_sum)
      type(ranger_tree), intent(in) :: tree
      real(dp), intent(in) :: x(:,:), case_weight(:)
      integer, intent(in) :: y(:), ncat(:), inbag(:)
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      real(dp), intent(inout) :: sum_imp(:), sumsq_imp(:)
      real(dp), intent(inout), optional :: local_sum(:,:)
      real(dp), allocatable :: xperm(:,:)
      integer, allocatable :: pred(:), base(:), idx(:), permutation(:)
      integer :: m, k, i, noob
      real(dp) :: base_err, perm_err, delta

      idx = pack([(i, i = 1, size(y))], [(is_evaluation_oob(i, inbag, case_weight, options%holdout), i = 1, size(y))])
      noob = size(idx)
      if (noob == 0) return
      allocate(xperm(size(x, 1), size(x, 2)), pred(size(y)), base(size(y)), permutation(noob))
      permutation = idx
      call predict_classification_tree(tree, x, ncat, prediction=base)
      base_err = real(count(base(idx) /= y(idx)), dp) / real(noob, dp)
      do m = 1, size(sum_imp)
         if (.not. tree_uses_variable(tree, m)) cycle
         perm_err = 0.0_dp
         do k = 1, max(1, options%n_perm)
            xperm = x
            call shuffle_int(rng, permutation)
            xperm(idx, m) = x(permutation, m)
            call predict_classification_tree(tree, xperm, ncat, prediction=pred)
            perm_err = perm_err + real(count(pred(idx) /= y(idx)), dp) / real(noob, dp)
            if (present(local_sum)) then
               do i = 1, noob
                  local_sum(m, idx(i)) = local_sum(m, idx(i)) + &
                     real(merge(1, 0, pred(idx(i)) /= y(idx(i))) - merge(1, 0, base(idx(i)) /= y(idx(i))), dp) / &
                     real(max(1, options%n_perm), dp)
               end do
            end if
         end do
         perm_err = perm_err / real(max(1, options%n_perm), dp)
         delta = perm_err - base_err
         sum_imp(m) = sum_imp(m) + delta
         sumsq_imp(m) = sumsq_imp(m) + delta * delta
      end do
   end subroutine classification_permutation_importance

   subroutine probability_permutation_importance(tree, x, y, ncat, inbag, case_weight, options, rng, &
      sum_imp, sumsq_imp, local_sum)
      type(ranger_tree), intent(in) :: tree
      real(dp), intent(in) :: x(:,:), case_weight(:)
      integer, intent(in) :: y(:), ncat(:), inbag(:)
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      real(dp), intent(inout) :: sum_imp(:), sumsq_imp(:)
      real(dp), intent(inout), optional :: local_sum(:,:)
      real(dp), allocatable :: xperm(:,:), pred(:,:), base(:,:)
      integer, allocatable :: idx(:), permutation(:)
      integer :: m, k, i, noob
      real(dp) :: base_err, perm_err, delta, base_case, perm_case

      idx = pack([(i, i = 1, size(y))], [(is_evaluation_oob(i, inbag, case_weight, options%holdout), i = 1, size(y))])
      noob = size(idx)
      if (noob == 0) return
      allocate(xperm(size(x, 1), size(x, 2)), pred(size(y), tree%nclass), base(size(y), tree%nclass), &
         permutation(noob))
      permutation = idx
      call predict_classification_tree(tree, x, ncat, probability=base)
      base_err = 0.0_dp
      do i = 1, noob
         base_err = base_err + (1.0_dp - base(idx(i), y(idx(i)))) ** 2
      end do
      base_err = base_err / real(noob, dp)
      do m = 1, size(sum_imp)
         if (.not. tree_uses_variable(tree, m)) cycle
         perm_err = 0.0_dp
         do k = 1, max(1, options%n_perm)
            xperm = x
            call shuffle_int(rng, permutation)
            xperm(idx, m) = x(permutation, m)
            call predict_classification_tree(tree, xperm, ncat, probability=pred)
            do i = 1, noob
               perm_case = (1.0_dp - pred(idx(i), y(idx(i)))) ** 2
               perm_err = perm_err + perm_case / real(max(1, options%n_perm), dp)
               if (present(local_sum)) then
                  base_case = (1.0_dp - base(idx(i), y(idx(i)))) ** 2
                  local_sum(m, idx(i)) = local_sum(m, idx(i)) + (perm_case - base_case) / &
                     real(max(1, options%n_perm), dp)
               end if
            end do
         end do
         perm_err = perm_err / real(noob, dp)
         delta = perm_err - base_err
         sum_imp(m) = sum_imp(m) + delta
         sumsq_imp(m) = sumsq_imp(m) + delta * delta
      end do
   end subroutine probability_permutation_importance

   subroutine finish_importance(sum_imp, sumsq_imp, ntree, scale, importance, importance_sd)
      real(dp), intent(in) :: sum_imp(:), sumsq_imp(:)
      integer, intent(in) :: ntree
      logical, intent(in) :: scale
      real(dp), intent(out) :: importance(:), importance_sd(:)
      real(dp) :: variance
      integer :: j
      importance = sum_imp / real(ntree, dp)
      importance_sd = 0.0_dp
      if (ntree > 0) then
         do j = 1, size(sum_imp)
            variance = max(0.0_dp, sumsq_imp(j) / real(ntree, dp) - importance(j) ** 2)
            importance_sd(j) = sqrt(variance / real(ntree, dp))
         end do
      end if
      if (scale) then
         do j = 1, size(importance)
            if (importance_sd(j) > 0.0_dp) importance(j) = importance(j) / importance_sd(j)
         end do
      end if
   end subroutine finish_importance

end module ranger_classification
