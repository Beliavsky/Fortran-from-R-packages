! SPDX-License-Identifier: GPL-3.0-only
! Modern Fortran translation of ranger 0.18.0 computational code.
! Upstream C++ core: Copyright (c) 2014-2018 Marvin N. Wright; MIT licensed.
! The upstream R package is GPL-3. See NOTICE.md for full provenance.
module ranger_tree_common
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp, i64
   use ranger_rng, only : ranger_rng_state, sample_indices
   use ranger_types, only : ranger_options, ranger_tree, RANGER_TERMINAL, RANGER_INTERIOR
   use ranger_types, only : RANGER_IMPORTANCE_IMPURITY_CORRECTED
   implicit none
   private

   real(dp), parameter :: split_eps = 64.0_dp * epsilon(1.0_dp)

   public :: initialize_tree, go_left_value, candidate_variables
   public :: collect_sorted_unique, partition_node, regularized_score
   public :: midpoint_safe, split_eps, update_best_split
   public :: base_predictor_count, unpermuted_var_id, is_shadow_variable
   public :: collect_present_categories, partition_mask_from_id, random_extratrees_factor_mask, tree_uses_variable

contains

   logical function tree_uses_variable(tree, var_id) result(used)
      type(ranger_tree), intent(in) :: tree
      integer, intent(in) :: var_id
      integer :: node

      used = .false.
      do node = 1, tree%n_nodes
         if (tree%status(node) /= RANGER_INTERIOR) cycle
         if (tree%split_var(node) == var_id) then
            used = .true.
            return
         end if
      end do
   end function tree_uses_variable

   subroutine initialize_tree(tree, maxnodes, maxcat, nvar, nclass, ntime)
      type(ranger_tree), intent(out) :: tree
      integer, intent(in) :: maxnodes, maxcat, nvar
      integer, intent(in), optional :: nclass, ntime
      integer :: nc, nt

      nc = 0
      nt = 0
      if (present(nclass)) nc = max(0, nclass)
      if (present(ntime)) nt = max(0, ntime)
      tree%n_nodes = 0
      tree%nclass = nc
      tree%ntime = nt
      tree%maxcat = max(1, maxcat)
      allocate(tree%left(maxnodes), tree%right(maxnodes), tree%split_var(maxnodes))
      allocate(tree%status(maxnodes), tree%node_class(maxnodes), tree%node_n(maxnodes))
      allocate(tree%split_value(maxnodes), tree%node_mean(maxnodes), tree%split_stat(maxnodes))
      allocate(tree%nan_go_right(maxnodes), tree%cat_left(tree%maxcat, maxnodes))
      allocate(tree%impurity_decrease(nvar))
      if (nc > 0) allocate(tree%class_prob(nc, maxnodes))
      if (nt > 0) allocate(tree%chf(nt, maxnodes))
      tree%left = 0
      tree%right = 0
      tree%split_var = 0
      tree%status = RANGER_TERMINAL
      tree%node_class = 0
      tree%node_n = 0
      tree%split_value = 0.0_dp
      tree%node_mean = 0.0_dp
      tree%split_stat = 0.0_dp
      tree%nan_go_right = .false.
      tree%cat_left = .false.
      tree%impurity_decrease = 0.0_dp
      if (nc > 0) tree%class_prob = 0.0_dp
      if (nt > 0) tree%chf = 0.0_dp
   end subroutine initialize_tree

   pure logical function go_left_value(value, ncat, cut, mask, nan_go_right) result(go_left)
      real(dp), intent(in) :: value, cut
      integer, intent(in) :: ncat
      logical, intent(in) :: mask(:), nan_go_right
      integer :: category

      if (ieee_is_nan(value)) then
         go_left = .not. nan_go_right
      else if (ncat <= 1) then
         go_left = value <= cut
      else
         category = nint(value)
         if (category < 1 .or. category > min(ncat, size(mask))) then
            go_left = .false.
         else
            go_left = mask(category)
         end if
      end if
   end function go_left_value

   pure integer function base_predictor_count(p, options) result(pbase)
      integer, intent(in) :: p
      type(ranger_options), intent(in) :: options
      pbase = p
      if (options%importance_mode == RANGER_IMPORTANCE_IMPURITY_CORRECTED .and. modulo(p, 2) == 0) pbase = p / 2
   end function base_predictor_count

   pure integer function unpermuted_var_id(var_id, p, options) result(base_id)
      integer, intent(in) :: var_id, p
      type(ranger_options), intent(in) :: options
      integer :: pbase
      pbase = base_predictor_count(p, options)
      base_id = var_id
      if (var_id > pbase) base_id = var_id - pbase
   end function unpermuted_var_id

   pure logical function is_shadow_variable(var_id, p, options) result(shadow)
      integer, intent(in) :: var_id, p
      type(ranger_options), intent(in) :: options
      shadow = var_id > base_predictor_count(p, options)
   end function is_shadow_variable


   subroutine collect_present_categories(values, indices, ncat, categories, ncategories)
      real(dp), intent(in) :: values(:)
      integer, intent(in) :: indices(:), ncat
      integer, allocatable, intent(out) :: categories(:)
      integer, intent(out) :: ncategories
      logical, allocatable :: seen(:)
      integer :: i, category

      allocate(seen(ncat), categories(ncat))
      seen = .false.
      ncategories = 0
      do i = 1, size(indices)
         if (ieee_is_nan(values(indices(i)))) cycle
         category = nint(values(indices(i)))
         if (category < 1 .or. category > ncat) cycle
         if (.not. seen(category)) then
            ncategories = ncategories + 1
            categories(ncategories) = category
            seen(category) = .true.
         end if
      end do
      if (ncategories > 1) call sort_integer_prefix(categories, ncategories)
   end subroutine collect_present_categories

   subroutine partition_mask_from_id(categories, ncategories, split_id, mask)
      integer, intent(in) :: categories(:), ncategories
      integer(i64), intent(in) :: split_id
      logical, intent(out) :: mask(:)
      integer :: j

      ! ranger stores the selected subset as the right child. The Fortran
      ! tree stores a logical left-child mask, so take the complement.
      mask = .true.
      do j = 1, ncategories
         if (btest(split_id, j - 1)) mask(categories(j)) = .false.
      end do
   end subroutine partition_mask_from_id

   subroutine random_extratrees_factor_mask(rng, values, indices, ncat, mask)
      type(ranger_rng_state), intent(inout) :: rng
      real(dp), intent(in) :: values(:)
      integer, intent(in) :: indices(:), ncat
      logical, intent(out) :: mask(:)
      logical, allocatable :: in_node(:), globally_observed(:)
      integer, allocatable :: present(:), absent(:)
      integer :: i, j, category, npresent, nabsent
      integer(i64) :: split_id, nchoices

      allocate(in_node(ncat), globally_observed(ncat), present(ncat), absent(ncat))
      in_node = .false.
      globally_observed = .false.

      do i = 1, size(values)
         if (ieee_is_nan(values(i))) cycle
         category = nint(values(i))
         if (category >= 1 .and. category <= ncat) globally_observed(category) = .true.
      end do
      do i = 1, size(indices)
         if (ieee_is_nan(values(indices(i)))) cycle
         category = nint(values(indices(i)))
         if (category >= 1 .and. category <= ncat) in_node(category) = .true.
      end do

      npresent = 0
      nabsent = 0
      do category = 1, ncat
         if (in_node(category)) then
            npresent = npresent + 1
            present(npresent) = category
         else if (globally_observed(category)) then
            nabsent = nabsent + 1
            absent(nabsent) = category
         end if
      end do

      ! Upstream splitID bits indicate the right child. Unselected and
      ! globally unused factor levels therefore go left.
      mask = .true.
      if (npresent > 1) then
         nchoices = shiftl(1_i64, npresent) - 2_i64
         split_id = rng%randint64(nchoices)
         do j = 1, npresent
            if (btest(split_id, j - 1)) mask(present(j)) = .false.
         end do
      end if
      if (nabsent > 1) then
         nchoices = shiftl(1_i64, nabsent)
         split_id = rng%randint64(nchoices) - 1_i64
         do j = 1, nabsent
            if (btest(split_id, j - 1)) mask(absent(j)) = .false.
         end do
      end if
   end subroutine random_extratrees_factor_mask

   subroutine candidate_variables(rng, p, mtry, options, vars, nvars)
      type(ranger_rng_state), intent(inout) :: rng
      integer, intent(in) :: p, mtry
      type(ranger_options), intent(in) :: options
      integer, allocatable, intent(out) :: vars(:)
      integer, intent(out) :: nvars
      integer, allocatable :: deterministic(:), available(:), draws(:), temp(:)
      logical, allocatable :: is_deterministic(:)
      real(dp), allocatable :: weights(:)
      integer :: i, j, k, ndet_original, ndet, navail, pbase, mrandom, rc
      logical :: corrected

      corrected = options%importance_mode == RANGER_IMPORTANCE_IMPURITY_CORRECTED .and. modulo(p, 2) == 0
      pbase = p
      if (corrected) pbase = p / 2
      allocate(deterministic(p), is_deterministic(p), temp(p))
      deterministic = 0
      is_deterministic = .false.
      ndet_original = 0
      if (allocated(options%always_split_variables)) then
         do i = 1, size(options%always_split_variables)
            j = options%always_split_variables(i)
            if (j < 1 .or. j > pbase) cycle
            if (is_deterministic(j)) cycle
            ndet_original = ndet_original + 1
            deterministic(ndet_original) = j
            is_deterministic(j) = .true.
         end do
      end if
      call sort_integer_prefix(deterministic, ndet_original)

      ndet = ndet_original
      if (corrected) then
         ! Preserve ranger 0.18.0 Forest::setAlwaysSplitVariables(): the
         ! shadow IDs use k + p, not original_varID + p. This is an
         ! upstream quirk and is observable when always-split variables
         ! are not the first variables.
         do k = 1, ndet_original
            j = pbase + k
            if (.not. is_deterministic(j)) then
               ndet = ndet + 1
               deterministic(ndet) = j
               is_deterministic(j) = .true.
            end if
         end do
         call sort_integer_prefix(deterministic, ndet)
      end if

      mrandom = max(0, mtry)
      allocate(draws(max(1, mrandom)))

      if (mrandom > 0) then
         if (allocated(options%split_select_weights)) then
            allocate(weights(p))
            weights = 0.0_dp
            if (size(options%split_select_weights) == pbase) then
               ! Match Forest::setSplitWeightVector() exactly: zero an
               ! original deterministic variable while ingesting the base
               ! weights, then copy the entire base half to the shadow half.
               ! The copy can therefore re-enable a shadow that is itself
               ! in deterministic_varIDs; createPossibleSplitVarSubset()
               ! does not remove deterministic IDs in the weighted branch.
               do i = 1, pbase
                  if (.not. is_deterministic(i)) weights(i) = options%split_select_weights(i)
               end do
               if (corrected) weights(pbase + 1:p) = weights(1:pbase)
            else
               weights = 1.0_dp
            end if
            call sample_indices(rng, p, mrandom, .false., draws(1:mrandom), weights=weights, status=rc)
            if (rc /= 0) mrandom = 0
         else
            allocate(available(p))
            navail = 0
            do i = 1, p
               if (.not. is_deterministic(i)) then
                  navail = navail + 1
                  available(navail) = i
               end if
            end do
            if (mrandom > navail) mrandom = navail
            if (mrandom > 0) then
               call sample_indices(rng, navail, mrandom, .false., draws(1:mrandom), status=rc)
               if (rc == 0) then
                  do i = 1, mrandom
                     draws(i) = available(draws(i))
                  end do
               else
                  mrandom = 0
               end if
            end if
         end if
      end if

      nvars = 0
      do i = 1, mrandom
         nvars = nvars + 1
         temp(nvars) = draws(i)
      end do
      do i = 1, ndet
         nvars = nvars + 1
         temp(nvars) = deterministic(i)
      end do
      allocate(vars(nvars))
      if (nvars > 0) vars = temp(1:nvars)
   end subroutine candidate_variables

   subroutine sort_integer_prefix(values, n)
      integer, intent(inout) :: values(:)
      integer, intent(in) :: n
      integer :: i, j, key

      do i = 2, n
         key = values(i)
         j = i - 1
         do while (j >= 1)
            if (values(j) <= key) exit
            values(j + 1) = values(j)
            j = j - 1
         end do
         values(j + 1) = key
      end do
   end subroutine sort_integer_prefix

   subroutine collect_sorted_unique(values, indices, unique_values, n_unique)
      real(dp), intent(in) :: values(:)
      integer, intent(in) :: indices(:)
      real(dp), allocatable, intent(out) :: unique_values(:)
      integer, intent(out) :: n_unique
      real(dp), allocatable :: work(:)
      integer :: i, j, nvalid
      real(dp) :: key

      nvalid = count(.not. ieee_is_nan(values(indices)))
      allocate(work(max(1, nvalid)))
      nvalid = 0
      do i = 1, size(indices)
         if (ieee_is_nan(values(indices(i)))) cycle
         nvalid = nvalid + 1
         work(nvalid) = values(indices(i))
      end do
      do i = 2, nvalid
         key = work(i)
         j = i - 1
         do while (j >= 1)
            if (work(j) <= key) exit
            work(j + 1) = work(j)
            j = j - 1
         end do
         work(j + 1) = key
      end do
      n_unique = 0
      do i = 1, nvalid
         if (i == 1) then
            n_unique = 1
            work(n_unique) = work(i)
         else if (abs(work(i) - work(n_unique)) > 0.0_dp) then
            n_unique = n_unique + 1
            work(n_unique) = work(i)
         end if
      end do
      allocate(unique_values(n_unique))
      if (n_unique > 0) unique_values = work(1:n_unique)
   end subroutine collect_sorted_unique

   subroutine partition_node(xcol, indices, weights, ncat, cut, mask, nan_go_right, left_idx, left_w, nleft, &
      right_idx, right_w, nright)
      real(dp), intent(in) :: xcol(:), weights(:), cut
      integer, intent(in) :: indices(:), ncat
      logical, intent(in) :: mask(:), nan_go_right
      integer, intent(out) :: left_idx(:), right_idx(:), nleft, nright
      real(dp), intent(out) :: left_w(:), right_w(:)
      integer :: i

      nleft = 0
      nright = 0
      do i = 1, size(indices)
         if (go_left_value(xcol(indices(i)), ncat, cut, mask, nan_go_right)) then
            nleft = nleft + 1
            left_idx(nleft) = indices(i)
            left_w(nleft) = weights(i)
         else
            nright = nright + 1
            right_idx(nright) = indices(i)
            right_w(nright) = weights(i)
         end if
      end do
   end subroutine partition_node

   real(dp) function regularized_score(score, var_id, depth, options, used_global, negative_score) result(value)
      real(dp), intent(in) :: score
      integer, intent(in) :: var_id, depth
      type(ranger_options), intent(in) :: options
      logical, intent(in) :: used_global(:)
      logical, intent(in), optional :: negative_score
      real(dp) :: factor
      logical :: negative

      value = score
      if (.not. allocated(options%regularization_factor)) return
      if (used_global(var_id)) return
      if (size(options%regularization_factor) == 1) then
         factor = options%regularization_factor(1)
      else if (size(options%regularization_factor) == size(used_global)) then
         factor = options%regularization_factor(var_id)
      else
         return
      end if
      if (factor >= 1.0_dp) return
      factor = max(factor, tiny(1.0_dp))
      if (options%regularization_usedepth) factor = factor ** real(depth + 1, dp)
      negative = .false.
      if (present(negative_score)) negative = negative_score
      if (negative) then
         value = score / factor
      else
         value = score * factor
      end if
   end function regularized_score

   pure real(dp) function midpoint_safe(a, b) result(mid)
      real(dp), intent(in) :: a, b
      mid = 0.5_dp * (a + b)
      if (abs(mid - b) <= 0.0_dp) mid = a
   end function midpoint_safe

   subroutine update_best_split(score, var_id, cut, mask, nan_right, best_score, best_var, best_cut, best_mask, &
      best_nan_right, found)
      real(dp), intent(in) :: score, cut
      integer, intent(in) :: var_id
      logical, intent(in) :: mask(:), nan_right
      real(dp), intent(inout) :: best_score, best_cut
      integer, intent(inout) :: best_var
      logical, intent(inout) :: best_mask(:), best_nan_right, found

      if (.not. found .or. score > best_score) then
         found = .true.
         best_score = score
         best_var = var_id
         best_cut = cut
         best_mask = .false.
         best_mask(1:min(size(mask), size(best_mask))) = mask(1:min(size(mask), size(best_mask)))
         best_nan_right = nan_right
      end if
   end subroutine update_best_split

end module ranger_tree_common
