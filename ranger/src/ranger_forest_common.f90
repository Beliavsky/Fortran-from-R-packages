! SPDX-License-Identifier: GPL-3.0-only
! Modern Fortran translation of ranger 0.18.0 computational code.
! Upstream C++ core: Copyright (c) 2014-2018 Marvin N. Wright; MIT licensed.
! The upstream R package is GPL-3. See NOTICE.md for full provenance.
module ranger_forest_common
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp, i64
   use ranger_rng, only : ranger_rng_state, sample_indices, shuffle_int
   use ranger_types, only : ranger_options, RANGER_IMPORTANCE_IMPURITY_CORRECTED, RANGER_SPLIT_MAXSTAT
   implicit none
   private

   public :: prepare_categories, make_case_weights, draw_inbag
   public :: draw_classwise_inbag, inbag_to_observations, default_sample_fraction
   public :: argmax_real, set_ok, set_error, validate_options, prepare_training_predictors
   public :: initialize_forest_rng, initialize_tree_rng, tree_options_for_index

contains

   subroutine initialize_forest_rng(rng, seed)
      type(ranger_rng_state), intent(inout) :: rng
      integer(i64), intent(in) :: seed
      integer(i64) :: clock_count

      if (seed == 0_i64) then
         call system_clock(count=clock_count)
         if (clock_count == 0_i64) clock_count = 1_i64
         call rng%seed(clock_count)
      else
         call rng%seed(seed)
      end if
   end subroutine initialize_forest_rng

   subroutine initialize_tree_rng(master_rng, tree_rng, seed, tree_index)
      type(ranger_rng_state), intent(inout) :: master_rng, tree_rng
      integer(i64), intent(in) :: seed
      integer, intent(in) :: tree_index
      integer(i64) :: tree_seed

      if (seed == 0_i64) then
         tree_seed = int(master_rng%uniform() * 4294967296.0_dp, i64)
      else
         tree_seed = modulo(seed * int(tree_index, i64), 4294967296_i64)
      end if
      call tree_rng%seed(tree_seed)
   end subroutine initialize_tree_rng


   subroutine prepare_categories(x, ncat_in, na_learn, ncat, status)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in), optional :: ncat_in(:)
      logical, intent(in) :: na_learn
      integer, allocatable, intent(out) :: ncat(:)
      integer, intent(out) :: status
      integer :: i, j, code

      status = 0
      allocate(ncat(size(x, 2)))
      ncat = 1
      if (present(ncat_in)) then
         if (size(ncat_in) /= size(x, 2) .or. any(ncat_in < 1)) then
            status = 1
            return
         end if
         ncat = ncat_in
      end if
      do j = 1, size(ncat)
         if (ncat(j) <= 1) cycle
         do i = 1, size(x, 1)
            if (ieee_is_nan(x(i, j))) then
               if (.not. na_learn) then
                  status = 2
                  return
               end if
               cycle
            end if
            code = nint(x(i, j))
            if (x(i, j) < real(code, dp) .or. x(i, j) > real(code, dp)) then
               status = 3
               return
            end if
            if (code < 1 .or. code > ncat(j)) then
               status = 4
               return
            end if
         end do
      end do
   end subroutine prepare_categories

   subroutine validate_options(options, p, status, message)
      type(ranger_options), intent(in) :: options
      integer, intent(in) :: p
      integer, intent(out) :: status
      character(len=*), intent(out) :: message
      logical, allocatable :: deterministic(:)
      integer :: i, j, ndet, navailable

      status = 0
      message = ''
      if (options%num_trees <= 0) then
         status = 1
         message = 'num_trees must be positive'
      else if (options%mtry > p) then
         status = 2
         message = 'mtry cannot exceed the number of predictors'
      else if ((options%sample_fraction >= 0.0_dp .and. options%sample_fraction <= 0.0_dp) .or. &
         options%sample_fraction > 1.0_dp) then
         status = 3
         message = 'sample_fraction must be in (0,1] or negative for the ranger default'
      else if (options%num_random_splits <= 0) then
         status = 4
         message = 'num_random_splits must be positive'
      else if (options%n_perm <= 0) then
         status = 5
         message = 'n_perm must be positive'
      else if (options%alpha < 0.0_dp .or. options%alpha > 1.0_dp) then
         status = 9
         message = 'alpha must be between 0 and 1'
      else if (options%minprop < 0.0_dp .or. options%minprop > 0.5_dp) then
         status = 10
         message = 'minprop must be between 0 and 0.5'
      else if (options%respect_unordered_factors < 0 .or. options%respect_unordered_factors > 2) then
         status = 11
         message = 'respect_unordered_factors must be ignore, order, or partition'
      else if (options%split_rule == RANGER_SPLIT_MAXSTAT .and. allocated(options%regularization_factor)) then
         if (any(options%regularization_factor < 1.0_dp)) then
            status = 8
            message = 'regularization cannot be used with maxstat splitrule'
         end if
      end if
      if (status /= 0) return

      allocate(deterministic(p))
      deterministic = .false.
      ndet = 0
      if (allocated(options%always_split_variables)) then
         do i = 1, size(options%always_split_variables)
            j = options%always_split_variables(i)
            if (j < 1 .or. j > p) then
               status = 12
               message = 'always_split_variables contains an invalid predictor index'
               return
            end if
            if (.not. deterministic(j)) then
               deterministic(j) = .true.
               ndet = ndet + 1
            end if
         end do
         if (ndet + options%mtry > p) then
            status = 13
            message = 'always-split variables plus mtry cannot exceed the number of predictors'
            return
         end if
      end if

      if (allocated(options%split_select_weights) .and. allocated(options%split_select_weights_by_tree)) then
         status = 15
         message = 'use either split_select_weights or split_select_weights_by_tree, not both'
         return
      end if
      if (allocated(options%split_select_weights)) then
         call validate_split_weights(options%split_select_weights, deterministic, p, options%mtry, &
            options%importance_mode, navailable)
         if (navailable < 0) then
            status = 6
            message = 'split_select_weights must be in [0,1] and have length nvar'
            return
         else if (navailable < options%mtry) then
            status = 14
            message = 'too many zero split-select weights for the requested mtry'
            return
         end if
      end if
      if (allocated(options%split_select_weights_by_tree)) then
         if (size(options%split_select_weights_by_tree, 1) /= options%num_trees .or. &
            size(options%split_select_weights_by_tree, 2) /= p) then
            status = 16
            message = 'split_select_weights_by_tree must have shape (num_trees,nvar)'
            return
         end if
         do j = 1, options%num_trees
            call validate_split_weights(options%split_select_weights_by_tree(j, :), deterministic, p, options%mtry, &
               options%importance_mode, navailable)
            if (navailable < 0) then
               status = 17
               message = 'per-tree split select weights must be in [0,1]'
               return
            else if (navailable < options%mtry) then
               status = 18
               message = 'a per-tree split weight row has too few positive selectable variables'
               return
            end if
         end do
      end if
      if (allocated(options%regularization_factor)) then
         if ((size(options%regularization_factor) /= 1 .and. size(options%regularization_factor) /= p) .or. &
            any(options%regularization_factor <= 0.0_dp) .or. any(options%regularization_factor > 1.0_dp)) then
            status = 7
            message = 'regularization_factor must be in (0,1] and have length 1 or nvar'
         end if
      end if
   end subroutine validate_options

   subroutine validate_split_weights(weights, deterministic, p, mtry, importance_mode, navailable)
      real(dp), intent(in) :: weights(:)
      logical, intent(in) :: deterministic(:)
      integer, intent(in) :: p, mtry, importance_mode
      integer, intent(out) :: navailable
      integer :: i

      navailable = -1
      if (size(weights) /= p .or. any(weights < 0.0_dp) .or. any(weights > 1.0_dp)) return
      navailable = 0
      do i = 1, p
         if (.not. deterministic(i) .and. weights(i) > 0.0_dp) navailable = navailable + 1
      end do
      if (importance_mode == RANGER_IMPORTANCE_IMPURITY_CORRECTED) then
         navailable = 2 * p - (p - navailable)
      end if
      if (mtry <= 0) navailable = max(navailable, 0)
   end subroutine validate_split_weights

   subroutine tree_options_for_index(options, tree_index, tree_options)
      type(ranger_options), intent(in) :: options
      integer, intent(in) :: tree_index
      type(ranger_options), intent(out) :: tree_options

      tree_options = options
      if (allocated(options%split_select_weights_by_tree)) then
         if (allocated(tree_options%split_select_weights)) deallocate(tree_options%split_select_weights)
         tree_options%split_select_weights = options%split_select_weights_by_tree(tree_index, :)
         deallocate(tree_options%split_select_weights_by_tree)
      end if
   end subroutine tree_options_for_index

   subroutine prepare_training_predictors(rng, x, ncat, options, x_train, ncat_train)
      type(ranger_rng_state), intent(in) :: rng
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: ncat(:)
      type(ranger_options), intent(in) :: options
      real(dp), allocatable, intent(out) :: x_train(:,:)
      integer, allocatable, intent(out) :: ncat_train(:)
      integer, allocatable :: permutation(:)
      type(ranger_rng_state) :: shadow_rng
      integer :: n, p, i

      n = size(x, 1)
      p = size(x, 2)
      if (options%importance_mode /= RANGER_IMPORTANCE_IMPURITY_CORRECTED) then
         allocate(x_train(n, p), ncat_train(p))
         x_train = x
         ncat_train = ncat
         return
      end if

      ! Upstream Data::permuteSampleIDs() creates one forest-wide shadow copy
      ! by applying the same random row permutation to every predictor.
      allocate(x_train(n, 2 * p), ncat_train(2 * p), permutation(n))
      x_train(:, 1:p) = x
      ncat_train(1:p) = ncat
      permutation = [(i, i = 1, n)]
      shadow_rng = rng
      call shuffle_int(shadow_rng, permutation)
      do i = 1, p
         x_train(:, p + i) = x(permutation, i)
         ncat_train(p + i) = ncat(i)
      end do
   end subroutine prepare_training_predictors

   real(dp) function default_sample_fraction(options) result(frac)
      type(ranger_options), intent(in) :: options
      frac = options%sample_fraction
      if (frac <= 0.0_dp) then
         if (options%replace) then
            frac = 1.0_dp
         else
            frac = 0.632_dp
         end if
      end if
   end function default_sample_fraction

   subroutine make_case_weights(n, case_weights, weights, status, active)
      integer, intent(in) :: n
      real(dp), intent(in), optional :: case_weights(:)
      real(dp), allocatable, intent(out) :: weights(:)
      integer, intent(out) :: status
      logical, intent(out), optional :: active

      status = 0
      if (present(active)) active = .false.
      allocate(weights(n))
      weights = 1.0_dp
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp)) then
            status = 1
            return
         end if
         if (sum(case_weights) <= 0.0_dp) then
            status = 2
            return
         end if
         ! ranger.R treats NULL or a constant vector as unweighted.
         if (all(abs(case_weights - case_weights(1)) <= 0.0_dp)) return
         weights = case_weights
         if (present(active)) active = .true.
      end if
   end subroutine make_case_weights

   subroutine draw_inbag(rng, options, case_weight, use_case_weights, inbag, status)
      type(ranger_rng_state), intent(inout) :: rng
      type(ranger_options), intent(in) :: options
      real(dp), intent(in) :: case_weight(:)
      logical, intent(in) :: use_case_weights
      integer, intent(out) :: inbag(:)
      integer, intent(out) :: status
      integer, allocatable :: sampled(:), pool(:)
      type(ranger_rng_state) :: shuffle_rng
      real(dp) :: frac
      integer :: sample_size, i

      status = 0
      inbag = 0
      frac = default_sample_fraction(options)
      if (options%holdout .and. use_case_weights) then
         sample_size = int(frac * real(count(case_weight > 0.0_dp), dp))
      else
         sample_size = int(frac * real(size(case_weight), dp))
      end if
      if (sample_size < 1) then
         status = 5
         return
      end if

      if (use_case_weights) then
         if (.not. options%replace .and. sample_size > count(case_weight > 0.0_dp)) then
            status = 6
            return
         end if
         allocate(sampled(sample_size))
         call sample_indices(rng, size(case_weight), sample_size, options%replace, sampled, weights=case_weight, status=status)
         if (status /= 0) return
         do i = 1, sample_size
            inbag(sampled(i)) = inbag(sampled(i)) + 1
         end do
      else if (options%replace) then
         allocate(sampled(sample_size))
         call sample_indices(rng, size(case_weight), sample_size, .true., sampled, status=status)
         if (status /= 0) return
         do i = 1, sample_size
            inbag(sampled(i)) = inbag(sampled(i)) + 1
         end do
      else
         ! ranger's shuffleAndSplit() takes the tree RNG by value.  The
         ! shuffle therefore selects the sample but does not advance the
         ! RNG subsequently used for split-variable selection.
         allocate(pool(size(case_weight)))
         pool = [(i, i = 1, size(case_weight))]
         shuffle_rng = rng
         call shuffle_int(shuffle_rng, pool)
         do i = 1, sample_size
            inbag(pool(i)) = 1
         end do
      end if
   end subroutine draw_inbag

   subroutine draw_classwise_inbag(rng, y, nclass, options, class_fraction, inbag, status)
      type(ranger_rng_state), intent(inout) :: rng
      integer, intent(in) :: y(:), nclass
      type(ranger_options), intent(in) :: options
      real(dp), intent(in) :: class_fraction(:)
      integer, intent(out) :: inbag(:)
      integer, intent(out) :: status
      integer, allocatable :: pool(:), sampled(:)
      type(ranger_rng_state) :: shuffle_rng
      integer :: k, i, nclass_obs, sample_size

      status = 0
      inbag = 0
      if (size(class_fraction) /= nclass) then
         status = 1
         return
      end if
      if (any(class_fraction < 0.0_dp) .or. any(class_fraction > 1.0_dp) .or. sum(class_fraction) <= 0.0_dp) then
         status = 2
         return
      end if
      do k = 1, nclass
         pool = pack([(i, i = 1, size(y))], y == k)
         nclass_obs = size(pool)
         if (nclass_obs == 0) cycle
         sample_size = nint(class_fraction(k) * real(size(y), dp))
         if (sample_size == 0) then
            deallocate(pool)
            cycle
         end if
         if (.not. options%replace .and. sample_size > nclass_obs) then
            status = 3
            return
         end if
         if (options%replace) then
            allocate(sampled(sample_size))
            call sample_indices(rng, nclass_obs, sample_size, .true., sampled, status=status)
            if (status /= 0) return
            do i = 1, sample_size
               inbag(pool(sampled(i))) = inbag(pool(sampled(i))) + 1
            end do
            deallocate(sampled)
         else
            ! shuffleAndSplitAppend() also takes the RNG by value. Each
            ! class shuffle starts from the same tree RNG state and leaves
            ! the tree RNG unchanged.
            shuffle_rng = rng
            call shuffle_int(shuffle_rng, pool)
            do i = 1, sample_size
               inbag(pool(i)) = 1
            end do
         end if
         deallocate(pool)
      end do
   end subroutine draw_classwise_inbag

   subroutine inbag_to_observations(inbag, obs_index, obs_weight)
      integer, intent(in) :: inbag(:)
      integer, allocatable, intent(out) :: obs_index(:)
      real(dp), allocatable, intent(out) :: obs_weight(:)
      integer :: i, j, k, nsampled

      ! Upstream ranger keeps bootstrap duplicates as distinct entries in
      ! Tree::sampleIDs.  Do the same here instead of compressing them into
      ! integer weights: sample counts affect min.node.size, maxstat ranks,
      ! and other node-size semantics even when impurity algebra is weighted.
      nsampled = sum(inbag)
      allocate(obs_index(nsampled), obs_weight(nsampled))
      j = 0
      do i = 1, size(inbag)
         do k = 1, inbag(i)
            j = j + 1
            obs_index(j) = i
            obs_weight(j) = 1.0_dp
         end do
      end do
   end subroutine inbag_to_observations

   integer function argmax_real(x) result(index)
      real(dp), intent(in) :: x(:)
      integer :: i
      index = 1
      do i = 2, size(x)
         if (x(i) > x(index)) index = i
      end do
   end function argmax_real

   subroutine set_ok(status, message)
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      if (present(status)) status = 0
      if (present(message)) message = ''
   end subroutine set_ok

   subroutine set_error(code, text, status, message)
      integer, intent(in) :: code
      character(len=*), intent(in) :: text
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      if (present(status)) status = code
      if (present(message)) message = text
   end subroutine set_error

end module ranger_forest_common
