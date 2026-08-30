! SPDX-License-Identifier: GPL-3.0-only
! Modern Fortran translation of ranger 0.18.0 computational code.
! Upstream C++ core: Copyright (c) 2014-2018 Marvin N. Wright; MIT licensed.
! The upstream R package is GPL-3. See NOTICE.md for full provenance.
module ranger_survival
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use ranger_rng, only : ranger_rng_state, shuffle_int, shuffle_real
   use ranger_types, only : ranger_options, ranger_survival_forest, ranger_tree
   use ranger_types, only : RANGER_IMPORTANCE_IMPURITY, RANGER_IMPORTANCE_IMPURITY_CORRECTED
   use ranger_types, only : RANGER_IMPORTANCE_PERMUTATION, RANGER_IMPORTANCE_PERMUTATION_CASEWISE
   use ranger_types, only : RANGER_SPLIT_STANDARD, RANGER_SPLIT_AUC, RANGER_SPLIT_AUC_IGNORE_TIES
   use ranger_types, only : RANGER_SPLIT_MAXSTAT, RANGER_SPLIT_EXTRATREES, RANGER_UNORDERED_IGNORE
   use ranger_types, only : RANGER_UNORDERED_ORDER, RANGER_UNORDERED_PARTITION
   use ranger_forest_common, only : prepare_categories, make_case_weights, draw_inbag, inbag_to_observations
   use ranger_forest_common, only : set_ok, set_error, validate_options, prepare_training_predictors
   use ranger_forest_common, only : initialize_forest_rng, initialize_tree_rng, tree_options_for_index
   use ranger_tree_survival, only : build_survival_tree, predict_survival_tree, concordance_index, concordance_casewise
   use ranger_tree_common, only : tree_uses_variable
   use ranger_factor_order, only : order_factors_survival, identity_factor_map, apply_factor_map
   implicit none
   private

   public :: fit_ranger_survival, predict_ranger_survival

   interface fit_ranger_survival
      module procedure fit_ranger_survival_vector
      module procedure fit_ranger_survival_count
      module procedure fit_ranger_survival_count_real
   end interface fit_ranger_survival

contains

   subroutine fit_ranger_survival_vector(x, time, event, forest, ncat, options, case_weights, time_interest, &
      preset_inbag, status, message)
      real(dp), intent(in) :: x(:,:), time(:)
      integer, intent(in) :: event(:)
      type(ranger_survival_forest), intent(out) :: forest
      integer, intent(in), optional :: ncat(:), preset_inbag(:,:)
      type(ranger_options), intent(in), optional :: options
      real(dp), intent(in), optional :: case_weights(:), time_interest(:)
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      type(ranger_options) :: opt, tree_opt
      type(ranger_rng_state) :: rng, tree_rng
      integer, allocatable :: cats(:), cats_train(:), inbag(:), obs_index(:), leaf(:), category_map(:,:)
      real(dp), allocatable :: obs_weight(:), case_weight(:), grid(:), tree_chf(:,:), imp_sum(:), imp_sumsq(:)
      real(dp), allocatable :: x_train(:,:), x_ordered(:,:)
      real(dp), allocatable :: local_sum(:,:)
      logical, allocatable :: used_global(:)
      logical :: use_case_weights
      integer :: n, p, t, i, rc
      character(len=200) :: message_local

      call set_ok(status, message)
      n = size(x, 1)
      p = size(x, 2)
      if (n <= 1 .or. p <= 0 .or. size(time) /= n .or. size(event) /= n) then
         call set_error(1, 'x, time, and event have incompatible or empty dimensions', status, message)
         return
      end if
      if (any((.not. ieee_is_finite(x)) .and. (.not. ieee_is_nan(x))) .or. any(.not. ieee_is_finite(time))) then
         call set_error(2, 'x contains infinite values or time contains non-finite values', status, message)
         return
      end if
      if (any(time < 0.0_dp) .or. any(event < 0) .or. any(event > 1) .or. count(event == 1) == 0) then
         call set_error(3, 'time must be nonnegative and event must be 0/1 with at least one event', status, message)
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
      if (opt%split_rule /= RANGER_SPLIT_STANDARD .and. opt%split_rule /= RANGER_SPLIT_AUC .and. &
         opt%split_rule /= RANGER_SPLIT_AUC_IGNORE_TIES .and. opt%split_rule /= RANGER_SPLIT_MAXSTAT .and. &
         opt%split_rule /= RANGER_SPLIT_EXTRATREES) then
         call set_error(18, 'unsupported survival split rule', status, message)
         return
      end if
      call prepare_categories(x, ncat, opt%na_learn, cats, rc)
      if (rc /= 0) then
         call set_error(20 + rc, 'invalid categorical predictor coding or missing-value policy', status, message)
         return
      end if
      if (opt%respect_unordered_factors == RANGER_UNORDERED_PARTITION .and. &
         any(cats > digits(1.0_dp))) then
         call set_error(38, 'unordered partition factors are limited to 53 observed levels', status, message)
         return
      end if
      if (opt%respect_unordered_factors == RANGER_UNORDERED_PARTITION .and. any(cats > 1) .and. &
         (opt%split_rule == RANGER_SPLIT_MAXSTAT .or. opt%split_rule == RANGER_SPLIT_AUC .or. &
          opt%split_rule == RANGER_SPLIT_AUC_IGNORE_TIES)) then
         call set_error(29, 'unordered partition splitting is unavailable for maxstat and AUC rules', status, message)
         return
      end if
      if (opt%respect_unordered_factors == RANGER_UNORDERED_ORDER) then
         call order_factors_survival(x, time, event, cats, category_map, x_ordered)
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
      call make_time_grid(time, event, time_interest, grid, rc)
      if (rc /= 0) then
         call set_error(40 + rc, 'time_interest must be finite and nonempty', status, message)
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

      call initialize_survival_forest(forest, n, p, cats, grid, opt)
      forest%category_map = category_map
      allocate(inbag(n), leaf(n), tree_chf(n, size(grid)), used_global(p), imp_sum(p), imp_sumsq(p))
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
         call build_survival_tree(x_train, time, event, obs_index, obs_weight, cats_train, grid, tree_opt, tree_rng, &
            used_global, forest%trees(t))
         if (opt%oob_error) then
            call predict_survival_tree(forest%trees(t), x_ordered, cats, chf=tree_chf, terminal_node=leaf)
            do i = 1, n
               if (.not. is_evaluation_oob(i, inbag, case_weight, opt%holdout)) cycle
               forest%oob_chf(i, :) = forest%oob_chf(i, :) + tree_chf(i, :)
               forest%oob_count(i) = forest%oob_count(i) + 1
            end do
            call update_survival_error(time, event, forest%oob_chf, forest%oob_count, forest%prediction_error(t))
         end if
         if (opt%importance_mode == RANGER_IMPORTANCE_IMPURITY .or. &
            opt%importance_mode == RANGER_IMPORTANCE_IMPURITY_CORRECTED) then
            imp_sum = imp_sum + forest%trees(t)%impurity_decrease
         else if (opt%importance_mode == RANGER_IMPORTANCE_PERMUTATION .or. &
            opt%importance_mode == RANGER_IMPORTANCE_PERMUTATION_CASEWISE) then
            if (allocated(local_sum)) then
               call survival_permutation_importance(forest%trees(t), x_ordered, time, event, cats, inbag, case_weight, &
                  opt, tree_rng, imp_sum, imp_sumsq, local_sum)
            else
               call survival_permutation_importance(forest%trees(t), x_ordered, time, event, cats, inbag, case_weight, &
                  opt, tree_rng, imp_sum, imp_sumsq)
            end if
         end if
         deallocate(obs_index, obs_weight)
      end do
      if (opt%oob_error) then
         do i = 1, n
            if (forest%oob_count(i) > 0) forest%oob_chf(i, :) = forest%oob_chf(i, :) / real(forest%oob_count(i), dp)
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
   end subroutine fit_ranger_survival_vector

   subroutine fit_ranger_survival_count(x, time, event, forest, ncat, options, case_weights, time_interest, &
      preset_inbag, status, message)
      real(dp), intent(in) :: x(:,:), time(:)
      integer, intent(in) :: event(:)
      type(ranger_survival_forest), intent(out) :: forest
      integer, intent(in), optional :: ncat(:), preset_inbag(:,:)
      type(ranger_options), intent(in), optional :: options
      real(dp), intent(in), optional :: case_weights(:)
      integer, intent(in) :: time_interest
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      real(dp), allocatable :: grid(:)
      integer :: rc

      call make_time_grid_count(time, event, time_interest, grid, rc)
      if (rc /= 0) then
         call set_error(41, 'time_interest count must be a positive integer', status, message)
         return
      end if
      call fit_ranger_survival_vector(x, time, event, forest, ncat, options, case_weights, grid, &
         preset_inbag, status, message)
   end subroutine fit_ranger_survival_count

   subroutine fit_ranger_survival_count_real(x, time, event, forest, ncat, options, case_weights, time_interest, &
      preset_inbag, status, message)
      real(dp), intent(in) :: x(:,:), time(:)
      integer, intent(in) :: event(:)
      type(ranger_survival_forest), intent(out) :: forest
      integer, intent(in), optional :: ncat(:), preset_inbag(:,:)
      type(ranger_options), intent(in), optional :: options
      real(dp), intent(in), optional :: case_weights(:)
      real(dp), intent(in) :: time_interest
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      integer :: npoint
      real(dp) :: fractional

      if (.not. ieee_is_finite(time_interest) .or. time_interest < 1.0_dp .or. &
         time_interest > real(huge(npoint), dp)) then
         call set_error(41, 'time_interest count must be a positive integer', status, message)
         return
      end if
      fractional = modulo(time_interest, 1.0_dp)
      if (fractional > 0.0_dp) then
         call set_error(41, 'time_interest count must be a positive integer', status, message)
         return
      end if
      npoint = nint(time_interest)
      call fit_ranger_survival_count(x, time, event, forest, ncat, options, case_weights, npoint, &
         preset_inbag, status, message)
   end subroutine fit_ranger_survival_count_real

   subroutine predict_ranger_survival(forest, x, chf, survival, per_tree_chf, terminal_nodes)
      type(ranger_survival_forest), intent(in) :: forest
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: chf(:,:)
      real(dp), intent(out), optional :: survival(:,:), per_tree_chf(:,:,:)
      integer, intent(out), optional :: terminal_nodes(:,:)
      real(dp), allocatable :: tree_chf(:,:), x_mapped(:,:)
      integer, allocatable :: leaf(:)
      integer :: n, t

      n = size(x, 1)
      if (size(x, 2) /= forest%nvar .or. size(chf, 1) /= n .or. size(chf, 2) /= forest%ntime) &
         error stop 'predict_ranger_survival: incompatible dimensions'
      if (present(survival)) then
         if (size(survival, 1) /= n .or. size(survival, 2) /= forest%ntime) &
            error stop 'predict_ranger_survival: survival has wrong shape'
      end if
      if (present(per_tree_chf)) then
         if (size(per_tree_chf, 1) /= n .or. size(per_tree_chf, 2) /= forest%ntime .or. &
            size(per_tree_chf, 3) /= size(forest%trees)) error stop 'predict_ranger_survival: per_tree_chf has wrong shape'
      end if
      if (present(terminal_nodes)) then
         if (size(terminal_nodes, 1) /= n .or. size(terminal_nodes, 2) /= size(forest%trees)) &
            error stop 'predict_ranger_survival: terminal_nodes has wrong shape'
      end if
      allocate(tree_chf(n, forest%ntime), leaf(n))
      call apply_factor_map(x, forest%ncat, forest%category_map, x_mapped)
      chf = 0.0_dp
      do t = 1, size(forest%trees)
         call predict_survival_tree(forest%trees(t), x_mapped, forest%ncat, chf=tree_chf, terminal_node=leaf)
         chf = chf + tree_chf
         if (present(per_tree_chf)) per_tree_chf(:, :, t) = tree_chf
         if (present(terminal_nodes)) terminal_nodes(:, t) = leaf
      end do
      chf = chf / real(size(forest%trees), dp)
      if (present(survival)) survival = exp(-chf)
   end subroutine predict_ranger_survival

   subroutine make_time_grid(time, event, time_interest, grid, status)
      real(dp), intent(in) :: time(:)
      integer, intent(in) :: event(:)
      real(dp), intent(in), optional :: time_interest(:)
      real(dp), allocatable, intent(out) :: grid(:)
      integer, intent(out) :: status

      status = 0
      if (present(time_interest)) then
         if (size(time_interest) == 0 .or. any(.not. ieee_is_finite(time_interest))) then
            status = 1
            allocate(grid(0))
            return
         end if
         ! ranger.R applies sort(unique(time.interest)) before entering the
         ! C++ core.  Preserve that interface behavior for explicit grids.
         call sorted_unique_real(time_interest, grid)
      else
         call event_time_grid(time, event, grid)
      end if
   end subroutine make_time_grid

   subroutine make_time_grid_count(time, event, npoint, grid, status)
      real(dp), intent(in) :: time(:)
      integer, intent(in) :: event(:), npoint
      real(dp), allocatable, intent(out) :: grid(:)
      real(dp), allocatable :: all_times(:)
      real(dp) :: position
      integer :: i, idx, nuniq
      integer, intent(out) :: status

      status = 0
      if (npoint < 1) then
         status = 1
         allocate(grid(0))
         return
      end if
      call event_time_grid(time, event, all_times)
      nuniq = size(all_times)
      if (nuniq <= npoint) then
         allocate(grid(nuniq))
         grid = all_times
         return
      end if

      allocate(grid(npoint))
      if (npoint == 1) then
         grid(1) = all_times(1)
         return
      end if
      do i = 1, npoint
         position = 1.0_dp + real(i - 1, dp) * real(nuniq - 1, dp) / real(npoint - 1, dp)
         idx = round_half_to_even(position)
         idx = max(1, min(nuniq, idx))
         grid(i) = all_times(idx)
      end do
   end subroutine make_time_grid_count

   subroutine event_time_grid(time, event, grid)
      real(dp), intent(in) :: time(:)
      integer, intent(in) :: event(:)
      real(dp), allocatable, intent(out) :: grid(:)
      real(dp), allocatable :: work(:)

      work = pack(time, event == 1)
      call sorted_unique_real(work, grid)
   end subroutine event_time_grid

   subroutine sorted_unique_real(values, unique_values)
      real(dp), intent(in) :: values(:)
      real(dp), allocatable, intent(out) :: unique_values(:)
      real(dp), allocatable :: work(:)
      real(dp) :: key
      integer :: i, j, nuniq

      allocate(work(size(values)))
      work = values
      do i = 2, size(work)
         key = work(i)
         j = i - 1
         do while (j >= 1)
            if (work(j) <= key) exit
            work(j + 1) = work(j)
            j = j - 1
         end do
         work(j + 1) = key
      end do
      nuniq = 0
      do i = 1, size(work)
         if (nuniq == 0) then
            nuniq = 1
            work(nuniq) = work(i)
         else if (work(i) > work(nuniq)) then
            nuniq = nuniq + 1
            work(nuniq) = work(i)
         end if
      end do
      allocate(unique_values(nuniq))
      if (nuniq > 0) unique_values = work(1:nuniq)
   end subroutine sorted_unique_real

   integer function round_half_to_even(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: fraction
      integer :: lower

      lower = floor(x)
      fraction = x - real(lower, dp)
      if (fraction < 0.5_dp) then
         value = lower
      else if (fraction > 0.5_dp) then
         value = lower + 1
      else if (mod(lower, 2) == 0) then
         value = lower
      else
         value = lower + 1
      end if
   end function round_half_to_even

   subroutine initialize_survival_forest(forest, n, p, ncat, grid, options)
      type(ranger_survival_forest), intent(out) :: forest
      integer, intent(in) :: n, p, ncat(:)
      real(dp), intent(in) :: grid(:)
      type(ranger_options), intent(in) :: options
      forest%nvar = p
      forest%nobs = n
      forest%ntime = size(grid)
      allocate(forest%ncat(p), forest%trees(options%num_trees), forest%unique_timepoints(size(grid)))
      allocate(forest%category_map(max(1, maxval(ncat)), p))
      forest%category_map = 0
      allocate(forest%oob_chf(n, size(grid)), forest%oob_count(n), forest%prediction_error(options%num_trees))
      allocate(forest%variable_importance(p), forest%variable_importance_sd(p))
      forest%ncat = ncat
      forest%unique_timepoints = grid
      forest%oob_chf = 0.0_dp
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
   end subroutine initialize_survival_forest

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

   subroutine update_survival_error(time, event, chf_sum, oob_count, error)
      real(dp), intent(in) :: time(:), chf_sum(:,:)
      integer, intent(in) :: event(:), oob_count(:)
      real(dp), intent(out) :: error
      integer, allocatable :: idx(:)
      real(dp), allocatable :: risk(:)
      integer :: i, nuse

      idx = pack([(i, i = 1, size(time))], oob_count > 0)
      nuse = size(idx)
      if (nuse <= 1) then
         error = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      allocate(risk(nuse))
      do i = 1, nuse
         risk(i) = sum(chf_sum(idx(i), :)) / real(oob_count(idx(i)), dp)
      end do
      error = 1.0_dp - concordance_index(time(idx), event(idx), risk)
   end subroutine update_survival_error

   subroutine survival_permutation_importance(tree, x, time, event, ncat, inbag, case_weight, options, rng, &
      sum_imp, sumsq_imp, local_sum)
      type(ranger_tree), intent(in) :: tree
      real(dp), intent(in) :: x(:,:), time(:), case_weight(:)
      integer, intent(in) :: event(:), ncat(:), inbag(:)
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      real(dp), intent(inout) :: sum_imp(:), sumsq_imp(:)
      real(dp), intent(inout), optional :: local_sum(:,:)
      real(dp), allocatable :: xperm(:,:), pred(:,:), base(:,:), base_risk(:), perm_risk(:)
      real(dp), allocatable :: base_case(:), perm_case(:)
      integer, allocatable :: idx(:), permutation(:)
      integer :: m, k, i, noob
      real(dp) :: base_err, perm_err, delta

      idx = pack([(i, i = 1, size(time))], [(is_evaluation_oob(i, inbag, case_weight, options%holdout), &
         i = 1, size(time))])
      noob = size(idx)
      if (noob <= 1) return
      allocate(xperm(size(x, 1), size(x, 2)), pred(size(x, 1), tree%ntime), &
         base(size(x, 1), tree%ntime), base_risk(noob), perm_risk(noob), permutation(noob))
      if (present(local_sum)) allocate(base_case(noob), perm_case(noob))
      permutation = idx
      call predict_survival_tree(tree, x, ncat, chf=base)
      do i = 1, noob
         base_risk(i) = sum(base(idx(i), :))
      end do
      base_err = 1.0_dp - concordance_index(time(idx), event(idx), base_risk)
      if (present(local_sum)) call concordance_casewise(time(idx), event(idx), base_risk, base_case)
      do m = 1, size(sum_imp)
         if (.not. tree_uses_variable(tree, m)) cycle
         perm_err = 0.0_dp
         do k = 1, max(1, options%n_perm)
            xperm = x
            call shuffle_int(rng, permutation)
            xperm(idx, m) = x(permutation, m)
            call predict_survival_tree(tree, xperm, ncat, chf=pred)
            do i = 1, noob
               perm_risk(i) = sum(pred(idx(i), :))
            end do
            perm_err = perm_err + 1.0_dp - concordance_index(time(idx), event(idx), perm_risk)
            if (present(local_sum)) then
               call concordance_casewise(time(idx), event(idx), perm_risk, perm_case)
               do i = 1, noob
                  local_sum(m, idx(i)) = local_sum(m, idx(i)) + (perm_case(i) - base_case(i)) / &
                     real(max(1, options%n_perm), dp)
               end do
            end if
         end do
         perm_err = perm_err / real(max(1, options%n_perm), dp)
         delta = perm_err - base_err
         sum_imp(m) = sum_imp(m) + delta
         sumsq_imp(m) = sumsq_imp(m) + delta * delta
      end do
   end subroutine survival_permutation_importance

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

end module ranger_survival
