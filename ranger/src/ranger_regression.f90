! SPDX-License-Identifier: GPL-3.0-only
! Modern Fortran translation of ranger 0.18.0 computational code.
! Upstream C++ core: Copyright (c) 2014-2018 Marvin N. Wright; MIT licensed.
! The upstream R package is GPL-3. See NOTICE.md for full provenance.
module ranger_regression
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use ranger_rng, only : ranger_rng_state, shuffle_int, shuffle_real, sample_indices
   use ranger_types, only : ranger_options, ranger_regression_forest, ranger_tree
   use ranger_types, only : RANGER_IMPORTANCE_IMPURITY, RANGER_IMPORTANCE_IMPURITY_CORRECTED
   use ranger_types, only : RANGER_IMPORTANCE_PERMUTATION, RANGER_IMPORTANCE_PERMUTATION_CASEWISE
   use ranger_types, only : RANGER_SPLIT_STANDARD, RANGER_SPLIT_MAXSTAT, RANGER_SPLIT_EXTRATREES
   use ranger_types, only : RANGER_SPLIT_BETA, RANGER_SPLIT_POISSON, RANGER_UNORDERED_IGNORE
   use ranger_types, only : RANGER_UNORDERED_ORDER, RANGER_UNORDERED_PARTITION
   use ranger_forest_common, only : prepare_categories, make_case_weights, draw_inbag, inbag_to_observations
   use ranger_forest_common, only : set_ok, set_error, validate_options, prepare_training_predictors
   use ranger_forest_common, only : initialize_forest_rng, initialize_tree_rng, tree_options_for_index
   use ranger_tree_regression, only : build_regression_tree, predict_regression_tree
   use ranger_tree_common, only : tree_uses_variable
   use ranger_factor_order, only : order_factors_regression, identity_factor_map, apply_factor_map
   implicit none
   private

   public :: fit_ranger_regression, predict_ranger_regression, predict_ranger_quantiles
   public :: predict_ranger_quantiles_oob

contains

   subroutine fit_ranger_regression(x, y, forest, ncat, options, case_weights, preset_inbag, status, message)
      real(dp), intent(in) :: x(:,:), y(:)
      type(ranger_regression_forest), intent(out) :: forest
      integer, intent(in), optional :: ncat(:), preset_inbag(:,:)
      type(ranger_options), intent(in), optional :: options
      real(dp), intent(in), optional :: case_weights(:)
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      type(ranger_options) :: opt, tree_opt
      type(ranger_rng_state) :: rng, tree_rng
      integer, allocatable :: cats(:), cats_train(:), inbag(:), obs_index(:), leaf(:), random_order(:), category_map(:,:)
      real(dp), allocatable :: obs_weight(:), case_weight(:), tree_pred(:), imp_sum(:), imp_sumsq(:), local_sum(:,:)
      real(dp), allocatable :: x_train(:,:), x_ordered(:,:)
      logical, allocatable :: used_global(:)
      logical :: use_case_weights
      integer :: n, p, t, i, rc
      character(len=200) :: message_local

      call set_ok(status, message)
      n = size(x, 1)
      p = size(x, 2)
      if (n <= 1 .or. p <= 0 .or. size(y) /= n) then
         call set_error(1, 'x and y have incompatible or empty dimensions', status, message)
         return
      end if
      if (any((.not. ieee_is_finite(x)) .and. (.not. ieee_is_nan(x))) .or. any(.not. ieee_is_finite(y))) then
         call set_error(2, 'x contains infinite values or y contains non-finite values', status, message)
         return
      end if
      opt = ranger_options()
      if (present(options)) opt = options
      if (opt%mtry <= 0) opt%mtry = max(1, int(sqrt(real(p, dp))))
      call validate_options(opt, p, rc, message_local)
      if (rc /= 0) then
         call set_error(10 + rc, trim(message_local), status, message)
         return
      end if
      if (opt%split_rule /= RANGER_SPLIT_STANDARD .and. opt%split_rule /= RANGER_SPLIT_MAXSTAT .and. &
         opt%split_rule /= RANGER_SPLIT_EXTRATREES .and. opt%split_rule /= RANGER_SPLIT_BETA .and. &
         opt%split_rule /= RANGER_SPLIT_POISSON) then
         call set_error(18, 'unsupported regression split rule', status, message)
         return
      end if
      if (opt%split_rule == RANGER_SPLIT_BETA) then
         if (any(y < 0.0_dp) .or. any(y > 1.0_dp)) then
            call set_error(20, 'beta splitrule requires response values between 0 and 1', status, message)
            return
         end if
      else if (opt%split_rule == RANGER_SPLIT_POISSON) then
         if (any(y < 0.0_dp) .or. sum(y) <= 0.0_dp) then
            call set_error(21, 'poisson splitrule requires y >= 0 with positive total response', status, message)
            return
         end if
      end if
      call prepare_categories(x, ncat, opt%na_learn, cats, rc)
      if (rc /= 0) then
         call set_error(30 + rc, 'invalid categorical predictor coding or missing-value policy', status, message)
         return
      end if
      if (opt%respect_unordered_factors == RANGER_UNORDERED_PARTITION .and. &
         any(cats > digits(1.0_dp))) then
         call set_error(38, 'unordered partition factors are limited to 53 observed levels', status, message)
         return
      end if
      if (opt%respect_unordered_factors == RANGER_UNORDERED_PARTITION .and. any(cats > 1) .and. &
         (opt%split_rule == RANGER_SPLIT_MAXSTAT .or. opt%split_rule == RANGER_SPLIT_BETA .or. &
          opt%split_rule == RANGER_SPLIT_POISSON)) then
         call set_error(39, 'unordered partition splitting is unavailable for maxstat, beta, and poisson', status, message)
         return
      end if
      if (opt%respect_unordered_factors == RANGER_UNORDERED_ORDER) then
         call order_factors_regression(x, y, cats, category_map, x_ordered)
         opt%respect_unordered_factors = RANGER_UNORDERED_IGNORE
      else
         call identity_factor_map(cats, category_map)
         call apply_factor_map(x, cats, category_map, x_ordered)
      end if
      call make_case_weights(n, case_weights, case_weight, rc, use_case_weights)
      if (rc /= 0) then
         call set_error(40 + rc, 'case_weights must be nonnegative with positive total weight', status, message)
         return
      end if
      if (opt%holdout .and. .not. use_case_weights) then
         call set_error(43, 'nonconstant case_weights are required in holdout mode', status, message)
         return
      end if
      if (present(preset_inbag)) then
         if (size(preset_inbag, 1) /= n .or. size(preset_inbag, 2) /= opt%num_trees .or. any(preset_inbag < 0)) then
            call set_error(50, 'preset_inbag must have shape (nobs,num_trees) and be nonnegative', status, message)
            return
         end if
      end if
      if (present(preset_inbag) .and. use_case_weights) then
         call set_error(51, 'case_weights and preset_inbag cannot be combined', status, message)
         return
      end if

      call initialize_regression_forest(forest, n, p, cats, opt)
      forest%category_map = category_map
      if (allocated(forest%training_y)) forest%training_y = y
      allocate(inbag(n), tree_pred(n), leaf(n), random_order(n), used_global(p), imp_sum(p), imp_sumsq(p))
      used_global = .false.
      imp_sum = 0.0_dp
      imp_sumsq = 0.0_dp
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
         else if (present(preset_inbag)) then
            inbag = preset_inbag(:, t)
            rc = 0
         else
            call draw_inbag(tree_rng, opt, case_weight, .false., inbag, rc)
         end if
         if (rc /= 0 .or. count(inbag > 0) == 0) then
            call set_error(60, 'sampling failed or produced an empty in-bag sample', status, message)
            return
         end if
         if (allocated(forest%inbag)) forest%inbag(:, t) = inbag
         call inbag_to_observations(inbag, obs_index, obs_weight)
         if (present(preset_inbag) .and. .not. use_case_weights) call shuffle_int(tree_rng, obs_index)
         call tree_options_for_index(opt, t, tree_opt)
         call build_regression_tree(x_train, y, obs_index, obs_weight, cats_train, tree_opt, tree_rng, used_global, forest%trees(t))
         if (opt%oob_error) then
            call predict_regression_tree(forest%trees(t), x_ordered, cats, prediction=tree_pred, terminal_node=leaf)
         else if (opt%quantreg) then
            call predict_regression_tree(forest%trees(t), x_ordered, cats, terminal_node=leaf)
         end if
         if (allocated(forest%training_terminal)) forest%training_terminal(:, t) = leaf
         if (allocated(forest%random_node_value)) then
            random_order = [(i, i = 1, n)]
            call shuffle_int(rng, random_order)
            do i = 1, n
               if (.not. forest%random_node_has_value(leaf(random_order(i)), t)) then
                  forest%random_node_value(leaf(random_order(i)), t) = y(random_order(i))
                  forest%random_node_has_value(leaf(random_order(i)), t) = .true.
               end if
            end do
         end if
         if (opt%oob_error) then
            do i = 1, n
               if (.not. is_evaluation_oob(i, inbag, case_weight, opt%holdout)) cycle
               forest%oob_prediction(i) = forest%oob_prediction(i) + tree_pred(i)
               forest%oob_count(i) = forest%oob_count(i) + 1
            end do
            call update_regression_error(y, forest%oob_prediction, forest%oob_count, forest%prediction_error(t))
         end if
         if (opt%importance_mode == RANGER_IMPORTANCE_IMPURITY .or. &
            opt%importance_mode == RANGER_IMPORTANCE_IMPURITY_CORRECTED) then
            imp_sum = imp_sum + forest%trees(t)%impurity_decrease
         else if (opt%importance_mode == RANGER_IMPORTANCE_PERMUTATION .or. &
            opt%importance_mode == RANGER_IMPORTANCE_PERMUTATION_CASEWISE) then
            if (allocated(local_sum)) then
               call regression_permutation_importance(forest%trees(t), x_ordered, y, cats, inbag, case_weight, opt, tree_rng, &
                  imp_sum, imp_sumsq, local_sum)
            else
               call regression_permutation_importance(forest%trees(t), x_ordered, y, cats, inbag, case_weight, opt, tree_rng, &
                  imp_sum, imp_sumsq)
            end if
         end if
         deallocate(obs_index, obs_weight)
      end do

      if (opt%oob_error) then
         do i = 1, n
            if (forest%oob_count(i) > 0) then
               forest%oob_prediction(i) = forest%oob_prediction(i) / real(forest%oob_count(i), dp)
            else
               forest%oob_prediction(i) = ieee_value(0.0_dp, ieee_quiet_nan)
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
      if (opt%quantreg .and. opt%keep_inbag) then
         call prepare_oob_quantile_values(forest, rng, rc)
         if (rc /= 0) then
            call set_error(70 + rc, 'too few usable OOB trees for quantile regression', status, message)
            return
         end if
      end if
   end subroutine fit_ranger_regression

   subroutine predict_ranger_regression(forest, x, prediction, per_tree, terminal_nodes)
      type(ranger_regression_forest), intent(in) :: forest
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: prediction(:)
      real(dp), intent(out), optional :: per_tree(:,:)
      integer, intent(out), optional :: terminal_nodes(:,:)
      real(dp), allocatable :: tree_pred(:), x_mapped(:,:)
      integer, allocatable :: leaf(:)
      integer :: n, t

      n = size(x, 1)
      if (size(x, 2) /= forest%nvar .or. size(prediction) /= n) &
         error stop 'predict_ranger_regression: incompatible dimensions'
      if (present(per_tree)) then
         if (size(per_tree, 1) /= n .or. size(per_tree, 2) /= size(forest%trees)) &
            error stop 'predict_ranger_regression: per_tree has wrong shape'
      end if
      if (present(terminal_nodes)) then
         if (size(terminal_nodes, 1) /= n .or. size(terminal_nodes, 2) /= size(forest%trees)) &
            error stop 'predict_ranger_regression: terminal_nodes has wrong shape'
      end if
      allocate(tree_pred(n), leaf(n))
      call apply_factor_map(x, forest%ncat, forest%category_map, x_mapped)
      prediction = 0.0_dp
      do t = 1, size(forest%trees)
         call predict_regression_tree(forest%trees(t), x_mapped, forest%ncat, prediction=tree_pred, terminal_node=leaf)
         prediction = prediction + tree_pred
         if (present(per_tree)) per_tree(:, t) = tree_pred
         if (present(terminal_nodes)) terminal_nodes(:, t) = leaf
      end do
      prediction = prediction / real(size(forest%trees), dp)
   end subroutine predict_ranger_regression

   subroutine predict_ranger_quantiles(forest, x, probabilities, quantiles)
      type(ranger_regression_forest), intent(in) :: forest
      real(dp), intent(in) :: x(:,:), probabilities(:)
      real(dp), intent(out) :: quantiles(:,:)
      integer, allocatable :: leaf(:)
      real(dp), allocatable :: values(:), x_mapped(:,:)
      integer :: nnew, ntree, i, t, m

      if (.not. allocated(forest%random_node_value) .or. .not. allocated(forest%random_node_has_value)) &
         error stop 'predict_ranger_quantiles: forest was not fit with options%quantreg = .true.'
      nnew = size(x, 1)
      ntree = size(forest%trees)
      if (size(x, 2) /= forest%nvar .or. size(quantiles, 1) /= nnew .or. &
         size(quantiles, 2) /= size(probabilities)) error stop 'predict_ranger_quantiles: incompatible dimensions'
      if (any(probabilities < 0.0_dp) .or. any(probabilities > 1.0_dp)) &
         error stop 'predict_ranger_quantiles: probabilities must be in [0,1]'
      allocate(leaf(nnew), values(ntree))
      call apply_factor_map(x, forest%ncat, forest%category_map, x_mapped)
      do i = 1, nnew
         do t = 1, ntree
            call predict_regression_tree(forest%trees(t), x_mapped(i:i, :), forest%ncat, terminal_node=leaf(1:1))
            if (.not. forest%random_node_has_value(leaf(1), t)) &
               error stop 'predict_ranger_quantiles: missing prepared terminal-node value'
            values(t) = forest%random_node_value(leaf(1), t)
         end do
         call sort_values(values)
         do m = 1, size(probabilities)
            quantiles(i, m) = empirical_quantile(values, probabilities(m))
         end do
      end do
   end subroutine predict_ranger_quantiles

   subroutine predict_ranger_quantiles_oob(forest, probabilities, quantiles)
      type(ranger_regression_forest), intent(in) :: forest
      real(dp), intent(in) :: probabilities(:)
      real(dp), intent(out) :: quantiles(:,:)
      real(dp), allocatable :: values(:)
      integer :: i, m

      if (.not. allocated(forest%random_node_value_oob)) &
         error stop 'predict_ranger_quantiles_oob: OOB quantile values were not prepared'
      if (size(quantiles, 1) /= forest%nobs .or. size(quantiles, 2) /= size(probabilities)) &
         error stop 'predict_ranger_quantiles_oob: incompatible output dimensions'
      if (any(probabilities < 0.0_dp) .or. any(probabilities > 1.0_dp)) &
         error stop 'predict_ranger_quantiles_oob: probabilities must be in [0,1]'
      allocate(values(forest%quantile_oob_count))
      do i = 1, forest%nobs
         values = forest%random_node_value_oob(i, :)
         call sort_values(values)
         do m = 1, size(probabilities)
            quantiles(i, m) = empirical_quantile(values, probabilities(m))
         end do
      end do
   end subroutine predict_ranger_quantiles_oob

   subroutine prepare_oob_quantile_values(forest, rng, status)
      type(ranger_regression_forest), intent(inout) :: forest
      type(ranger_rng_state), intent(inout) :: rng
      integer, intent(out) :: status
      real(dp), allocatable :: candidates(:)
      integer, allocatable :: pool(:), selected(:)
      integer :: n, ntree, minoob, i, j, t, ncandidate, npool, rc, leaf

      status = 0
      if (.not. allocated(forest%inbag) .or. .not. allocated(forest%training_terminal) .or. &
         .not. allocated(forest%training_y)) then
         status = 1
         return
      end if
      n = forest%nobs
      ntree = size(forest%trees)
      minoob = ntree
      do i = 1, n
         minoob = min(minoob, count(forest%inbag(i, :) == 0))
      end do
      if (minoob < 10) then
         status = 2
         return
      end if
      forest%quantile_oob_count = minoob
      allocate(forest%random_node_value_oob(n, minoob), candidates(ntree), pool(n), selected(minoob))
      do i = 1, n
         ncandidate = 0
         do t = 1, ntree
            if (forest%inbag(i, t) > 0) cycle
            leaf = forest%training_terminal(i, t)
            npool = 0
            do j = 1, n
               if (j == i) cycle
               if (forest%training_terminal(j, t) == leaf) then
                  npool = npool + 1
                  pool(npool) = j
               end if
            end do
            if (npool == 0) cycle
            ncandidate = ncandidate + 1
            candidates(ncandidate) = forest%training_y(pool(rng%randint(npool)))
         end do
         if (ncandidate < minoob) then
            status = 3
            return
         end if
         call sample_indices(rng, ncandidate, minoob, .false., selected, status=rc)
         if (rc /= 0) then
            status = 4
            return
         end if
         forest%random_node_value_oob(i, :) = candidates(selected)
      end do
   end subroutine prepare_oob_quantile_values

   subroutine initialize_regression_forest(forest, n, p, ncat, options)
      type(ranger_regression_forest), intent(out) :: forest
      integer, intent(in) :: n, p, ncat(:)
      type(ranger_options), intent(in) :: options
      forest%nvar = p
      forest%nobs = n
      allocate(forest%ncat(p), forest%trees(options%num_trees), forest%oob_prediction(n), forest%oob_count(n))
      allocate(forest%category_map(max(1, maxval(ncat)), p))
      forest%category_map = 0
      allocate(forest%prediction_error(options%num_trees), forest%variable_importance(p), forest%variable_importance_sd(p))
      forest%ncat = ncat
      forest%oob_prediction = 0.0_dp
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
      if (options%quantreg) then
         allocate(forest%training_y(n), forest%training_terminal(n, options%num_trees))
         allocate(forest%random_node_value(2 * n + 1, options%num_trees))
         allocate(forest%random_node_has_value(2 * n + 1, options%num_trees))
         forest%training_terminal = 0
         forest%random_node_value = 0.0_dp
         forest%random_node_has_value = .false.
      end if
   end subroutine initialize_regression_forest

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

   subroutine update_regression_error(y, prediction_sum, oob_count, error)
      real(dp), intent(in) :: y(:), prediction_sum(:)
      integer, intent(in) :: oob_count(:)
      real(dp), intent(out) :: error
      integer :: i, nuse
      real(dp) :: pred, error_sum
      nuse = 0
      error = ieee_value(0.0_dp, ieee_quiet_nan)
      error_sum = 0.0_dp
      do i = 1, size(y)
         if (oob_count(i) <= 0) cycle
         pred = prediction_sum(i) / real(oob_count(i), dp)
         error_sum = error_sum + (pred - y(i)) ** 2
         nuse = nuse + 1
      end do
      if (nuse > 0) error = error_sum / real(nuse, dp)
   end subroutine update_regression_error

   subroutine regression_permutation_importance(tree, x, y, ncat, inbag, case_weight, options, rng, &
      sum_imp, sumsq_imp, local_sum)
      type(ranger_tree), intent(in) :: tree
      real(dp), intent(in) :: x(:,:), y(:), case_weight(:)
      integer, intent(in) :: ncat(:), inbag(:)
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      real(dp), intent(inout) :: sum_imp(:), sumsq_imp(:)
      real(dp), intent(inout), optional :: local_sum(:,:)
      real(dp), allocatable :: xperm(:,:), pred(:), base(:)
      integer, allocatable :: idx(:), permutation(:)
      integer :: m, k, i, noob
      real(dp) :: base_err, perm_err, delta, base_case, perm_case

      idx = pack([(i, i = 1, size(y))], [(is_evaluation_oob(i, inbag, case_weight, options%holdout), i = 1, size(y))])
      noob = size(idx)
      if (noob == 0) return
      allocate(xperm(size(x, 1), size(x, 2)), pred(size(y)), base(size(y)), permutation(noob))
      permutation = idx
      call predict_regression_tree(tree, x, ncat, prediction=base)
      base_err = sum((base(idx) - y(idx)) ** 2) / real(noob, dp)
      do m = 1, size(sum_imp)
         if (.not. tree_uses_variable(tree, m)) cycle
         perm_err = 0.0_dp
         do k = 1, max(1, options%n_perm)
            xperm = x
            call shuffle_int(rng, permutation)
            xperm(idx, m) = x(permutation, m)
            call predict_regression_tree(tree, xperm, ncat, prediction=pred)
            perm_err = perm_err + sum((pred(idx) - y(idx)) ** 2) / real(noob, dp)
            if (present(local_sum)) then
               do i = 1, noob
                  base_case = (base(idx(i)) - y(idx(i))) ** 2
                  perm_case = (pred(idx(i)) - y(idx(i))) ** 2
                  local_sum(m, idx(i)) = local_sum(m, idx(i)) + (perm_case - base_case) / &
                     real(max(1, options%n_perm), dp)
               end do
            end if
         end do
         perm_err = perm_err / real(max(1, options%n_perm), dp)
         delta = perm_err - base_err
         sum_imp(m) = sum_imp(m) + delta
         sumsq_imp(m) = sumsq_imp(m) + delta * delta
      end do
   end subroutine regression_permutation_importance

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

   subroutine sort_values(values)
      real(dp), intent(inout) :: values(:)
      integer :: i, j
      real(dp) :: value
      do i = 2, size(values)
         value = values(i)
         j = i - 1
         do while (j >= 1)
            if (values(j) <= value) exit
            values(j + 1) = values(j)
            j = j - 1
         end do
         values(j + 1) = value
      end do
   end subroutine sort_values

   real(dp) function empirical_quantile(values, probability) result(q)
      real(dp), intent(in) :: values(:), probability
      real(dp) :: h, frac
      integer :: lo, hi, n
      n = size(values)
      if (n == 1) then
         q = values(1)
         return
      end if
      h = 1.0_dp + real(n - 1, dp) * probability
      lo = max(1, min(n, int(floor(h))))
      hi = max(1, min(n, lo + 1))
      frac = h - real(lo, dp)
      q = (1.0_dp - frac) * values(lo) + frac * values(hi)
   end function empirical_quantile

end module ranger_regression
