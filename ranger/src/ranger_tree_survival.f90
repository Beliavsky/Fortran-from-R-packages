! SPDX-License-Identifier: GPL-3.0-only
! Modern Fortran translation of ranger 0.18.0 computational code.
! Upstream C++ core: Copyright (c) 2014-2018 Marvin N. Wright; MIT licensed.
! The upstream R package is GPL-3. See NOTICE.md for full provenance.
module ranger_tree_survival
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp, i64
   use ranger_rng, only : ranger_rng_state
   use ranger_types, only : ranger_options, ranger_tree, RANGER_INTERIOR
   use ranger_types, only : RANGER_SPLIT_STANDARD, RANGER_SPLIT_AUC, RANGER_SPLIT_AUC_IGNORE_TIES
   use ranger_types, only : RANGER_SPLIT_MAXSTAT, RANGER_SPLIT_EXTRATREES, RANGER_UNORDERED_PARTITION
   use ranger_tree_common, only : initialize_tree, go_left_value, candidate_variables, collect_sorted_unique
   use ranger_tree_common, only : midpoint_safe, regularized_score
   use ranger_tree_common, only : collect_present_categories, partition_mask_from_id
   use ranger_tree_common, only : random_extratrees_factor_mask
   use ranger_tree_common, only : base_predictor_count, unpermuted_var_id, is_shadow_variable
   use ranger_maxstat, only : logrank_scores, maxstat_best, maxstat_pvalue_lau92, maxstat_pvalue_lau94
   use ranger_maxstat, only : maxstat_pvalue_unadjusted, adjust_pvalues_bh
   implicit none
   private

   public :: build_survival_tree, predict_survival_tree, concordance_index, concordance_casewise

contains

   subroutine build_survival_tree(x, time, status, obs_index, obs_weight, ncat, timegrid, options, rng, &
      used_global, tree)
      real(dp), intent(in) :: x(:,:), time(:), obs_weight(:), timegrid(:)
      integer, intent(in) :: status(:), obs_index(:), ncat(:)
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      logical, intent(inout) :: used_global(:)
      type(ranger_tree), intent(out) :: tree
      integer :: maxnodes, maxcat, next_node, min_node, min_bucket, pbase

      pbase = base_predictor_count(size(x, 2), options)
      maxnodes = max(3, 2 * size(obs_index) + 1)
      maxcat = max(1, maxval(ncat))
      min_node = options%min_node_size
      if (min_node <= 0) min_node = 3
      min_bucket = options%min_bucket
      if (min_bucket <= 0) min_bucket = 3
      call initialize_tree(tree, maxnodes, maxcat, pbase, ntime=size(timegrid))
      next_node = 1
      call grow_node(obs_index, obs_weight, 1, 0)
      tree%n_nodes = next_node

   contains

      recursive subroutine grow_node(indices, weights, node, depth)
         integer, intent(in) :: indices(:), node, depth
         real(dp), intent(in) :: weights(:)
         integer, allocatable :: left_idx(:), right_idx(:)
         real(dp), allocatable :: left_w(:), right_w(:)
         logical, allocatable :: best_mask(:)
         integer :: best_var, base_var, nleft, nright, left_node, right_node
         real(dp) :: best_cut, best_score, importance_sign
         logical :: best_nan_right, found

         tree%node_n(node) = nint(sum(weights))
         call compute_node_chf(time, status, indices, weights, timegrid, tree%chf(:, node))
         tree%node_mean(node) = sum(tree%chf(:, node))

         if (size(indices) <= min_node) return
         if (options%max_depth > 0 .and. depth >= options%max_depth) return
         if (count_events(status, indices, weights) <= 0.0_dp) return
         if (next_node + 2 > size(tree%left)) return

         allocate(best_mask(maxcat))
         call best_survival_split(x, time, status, indices, weights, ncat, options, min_bucket, depth, rng, &
            used_global, best_var, best_cut, best_mask, best_nan_right, best_score, found)
         if (.not. found) return

         allocate(left_idx(size(indices)), right_idx(size(indices)))
         allocate(left_w(size(indices)), right_w(size(indices)))
         call split_indices(x(:, best_var), indices, weights, ncat(best_var), best_cut, best_mask, best_nan_right, &
            left_idx, left_w, nleft, right_idx, right_w, nright)
         if (nleft == 0 .or. nright == 0) return

         base_var = unpermuted_var_id(best_var, size(x, 2), options)
         importance_sign = 1.0_dp
         if (is_shadow_variable(best_var, size(x, 2), options)) importance_sign = -1.0_dp
         tree%status(node) = RANGER_INTERIOR
         tree%split_var(node) = base_var
         tree%split_value(node) = best_cut
         tree%split_stat(node) = best_score
         tree%nan_go_right(node) = best_nan_right
         if (ncat(best_var) > 1) tree%cat_left(1:ncat(best_var), node) = best_mask(1:ncat(best_var))
         tree%impurity_decrease(base_var) = tree%impurity_decrease(base_var) + importance_sign * best_score
         used_global(base_var) = .true.

         left_node = next_node + 1
         right_node = next_node + 2
         next_node = next_node + 2
         tree%left(node) = left_node
         tree%right(node) = right_node
         call grow_node(left_idx(1:nleft), left_w(1:nleft), left_node, depth + 1)
         call grow_node(right_idx(1:nright), right_w(1:nright), right_node, depth + 1)
      end subroutine grow_node

   end subroutine build_survival_tree

   subroutine compute_node_chf(time, status, indices, weights, grid, chf)
      real(dp), intent(in) :: time(:), weights(:), grid(:)
      integer, intent(in) :: status(:), indices(:)
      real(dp), intent(out) :: chf(size(grid))
      real(dp) :: at_risk, deaths, cumulative
      integer :: g, i

      cumulative = 0.0_dp
      chf = 0.0_dp
      do g = 1, size(grid)
         at_risk = 0.0_dp
         deaths = 0.0_dp
         do i = 1, size(indices)
            if (time(indices(i)) >= grid(g)) at_risk = at_risk + weights(i)
            if (status(indices(i)) == 1 .and. abs(time(indices(i)) - grid(g)) <= 0.0_dp) deaths = deaths + weights(i)
         end do
         if (at_risk > 0.0_dp .and. deaths > 0.0_dp) cumulative = cumulative + deaths / at_risk
         chf(g) = cumulative
      end do
   end subroutine compute_node_chf

   pure real(dp) function count_events(status, indices, weights) result(events)
      integer, intent(in) :: status(:), indices(:)
      real(dp), intent(in) :: weights(:)
      integer :: i
      events = 0.0_dp
      do i = 1, size(indices)
         if (status(indices(i)) == 1) events = events + weights(i)
      end do
   end function count_events

   subroutine best_survival_split(x, time, status, indices, weights, ncat, options, min_bucket, depth, rng, &
      used_global, best_var, best_cut, best_mask, best_nan_right, best_score, found)
      real(dp), intent(in) :: x(:,:), time(:), weights(:)
      integer, intent(in) :: status(:), indices(:), ncat(:), min_bucket, depth
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      logical, intent(in) :: used_global(:)
      integer, intent(out) :: best_var
      real(dp), intent(out) :: best_cut, best_score
      logical, intent(out) :: best_mask(:), best_nan_right, found
      integer, allocatable :: vars(:)
      logical, allocatable :: mask(:)
      integer :: i, v, base_v, nvars, mtry, pbase
      real(dp) :: score, cut, adjusted, pvalue, best_p
      logical :: nan_right, ok

      if (options%split_rule == RANGER_SPLIT_MAXSTAT) then
         call best_survival_split_maxstat(x, time, status, indices, ncat, options, rng, best_var, best_cut, &
            best_mask, best_nan_right, best_score, found)
         return
      end if

      pbase = base_predictor_count(size(x, 2), options)
      mtry = options%mtry
      if (mtry <= 0) mtry = max(1, int(sqrt(real(pbase, dp))))
      call candidate_variables(rng, size(x, 2), mtry, options, vars, nvars)
      allocate(mask(size(best_mask)))
      best_var = 0
      best_cut = 0.0_dp
      best_score = -huge(1.0_dp)
      best_mask = .false.
      best_nan_right = .false.
      best_p = 1.0_dp
      found = .false.
      do i = 1, nvars
         v = vars(i)
         mask = .false.
         call best_survival_variable(x(:, v), time, status, indices, weights, ncat(v), options, min_bucket, rng, &
            score, cut, mask, nan_right, pvalue, ok)
         if (.not. ok) cycle
         base_v = unpermuted_var_id(v, size(x, 2), options)
         adjusted = regularized_score(score, base_v, depth, options, used_global)
         if (.not. found .or. adjusted > best_score) then
            found = .true.
            best_var = v
            best_cut = cut
            best_mask = mask
            best_nan_right = nan_right
            best_score = adjusted
            best_p = pvalue
         end if
      end do
   end subroutine best_survival_split

   subroutine best_survival_split_maxstat(x, time, status, indices, ncat, options, rng, best_var, best_cut, &
      best_mask, best_nan_right, best_stat, found)
      real(dp), intent(in) :: x(:,:), time(:)
      integer, intent(in) :: status(:), indices(:), ncat(:)
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      integer, intent(out) :: best_var
      real(dp), intent(out) :: best_cut, best_stat
      logical, intent(out) :: best_mask(:), best_nan_right, found
      integer, allocatable :: vars(:), num_left(:), trimmed_left(:), map(:)
      real(dp), allocatable :: node_time(:), scores(:), xnode(:), pvalue(:), adjusted(:), stat(:), cut(:), pv(:), av(:)
      integer, allocatable :: node_status(:)
      logical, allocatable :: valid(:)
      integer :: nvars, mtry, pbase, i, v, nvalid, imin, c, j, ncuts
      real(dp) :: p92, p94
      logical :: ok

      pbase = base_predictor_count(size(x, 2), options)
      mtry = options%mtry
      if (mtry <= 0) mtry = max(1, int(sqrt(real(pbase, dp))))
      call candidate_variables(rng, size(x, 2), mtry, options, vars, nvars)
      allocate(node_time(size(indices)), node_status(size(indices)), scores(size(indices)), xnode(size(indices)))
      allocate(pvalue(nvars), adjusted(nvars), stat(nvars), cut(nvars), valid(nvars))
      node_time = time(indices)
      node_status = status(indices)
      call logrank_scores(node_time, node_status, scores)
      pvalue = 1.0_dp
      adjusted = 1.0_dp
      stat = -1.0_dp
      cut = 0.0_dp
      valid = .false.

      do i = 1, nvars
         v = vars(i)
         xnode = x(indices, v)
         call maxstat_best(scores, xnode, options%minprop, 1.0_dp - options%minprop, stat(i), cut(i), num_left, ok)
         if (.not. ok) cycle
         ncuts = max(0, size(num_left) - 1)
         if (ncuts == 1) then
            pvalue(i) = maxstat_pvalue_unadjusted(stat(i))
         else if (ncuts > 1) then
            allocate(trimmed_left(ncuts))
            trimmed_left = num_left(1:ncuts)
            p92 = maxstat_pvalue_lau92(stat(i), options%minprop, 1.0_dp - options%minprop)
            p94 = maxstat_pvalue_lau94(stat(i), options%minprop, 1.0_dp - options%minprop, size(indices), trimmed_left)
            pvalue(i) = min(p92, p94)
            deallocate(trimmed_left)
         else
            deallocate(num_left)
            cycle
         end if
         valid(i) = .true.
         deallocate(num_left)
      end do

      nvalid = count(valid)
      found = .false.
      best_var = 0
      best_cut = 0.0_dp
      best_stat = -huge(1.0_dp)
      best_mask = .false.
      best_nan_right = .false.
      if (nvalid == 0) return
      allocate(pv(nvalid), av(nvalid), map(nvalid))
      j = 0
      do i = 1, nvars
         if (.not. valid(i)) cycle
         j = j + 1
         pv(j) = pvalue(i)
         map(j) = i
      end do
      call adjust_pvalues_bh(pv, av)
      do j = 1, nvalid
         adjusted(map(j)) = av(j)
      end do

      imin = 0
      do i = 1, nvars
         if (.not. valid(i)) cycle
         if (imin == 0) then
            imin = i
         else if (pvalue(i) < pvalue(imin)) then
            imin = i
         end if
      end do
      if (imin == 0) return
      if (adjusted(imin) > options%alpha) return

      best_var = vars(imin)
      best_cut = cut(imin)
      best_stat = stat(imin)
      found = .true.
      if (ncat(best_var) > 1) then
         do c = 1, min(ncat(best_var), size(best_mask))
            best_mask(c) = real(c, dp) <= best_cut
         end do
      end if
   end subroutine best_survival_split_maxstat

   subroutine best_survival_variable(values, time, status, indices, weights, ncat, options, min_bucket, rng, &
      best_score, best_cut, best_mask, best_nan_right, best_pvalue, found)
      real(dp), intent(in) :: values(:), time(:), weights(:)
      integer, intent(in) :: status(:), indices(:), ncat, min_bucket
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      real(dp), intent(out) :: best_score, best_cut, best_pvalue
      logical, intent(out) :: best_mask(:), best_nan_right, found

      if (ncat <= 1 .or. options%respect_unordered_factors /= RANGER_UNORDERED_PARTITION) then
         call best_numeric_survival(values, time, status, indices, weights, options, min_bucket, rng, best_score, &
            best_cut, best_nan_right, best_pvalue, found)
         best_mask = .false.
      else
         call best_categorical_survival(values, time, status, indices, weights, ncat, options, min_bucket, rng, &
            best_score, best_mask, best_nan_right, best_pvalue, found)
         best_cut = 0.0_dp
      end if
   end subroutine best_survival_variable

   subroutine best_numeric_survival(values, time, status, indices, weights, options, min_bucket, rng, best_score, &
      best_cut, best_nan_right, best_pvalue, found)
      real(dp), intent(in) :: values(:), time(:), weights(:)
      integer, intent(in) :: status(:), indices(:), min_bucket
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      real(dp), intent(out) :: best_score, best_cut, best_pvalue
      logical, intent(out) :: best_nan_right, found
      real(dp), allocatable :: unique(:)
      integer :: nu, i, r, ntry
      real(dp) :: cut, lo, hi, score_l, score_r, p_l, p_r
      logical :: ok_l, ok_r

      call collect_sorted_unique(values, indices, unique, nu)
      best_score = -huge(1.0_dp)
      best_cut = 0.0_dp
      best_nan_right = .false.
      best_pvalue = 1.0_dp
      found = .false.
      if (nu < 2) return
      if (options%split_rule == RANGER_SPLIT_EXTRATREES) then
         ntry = max(1, options%num_random_splits)
         lo = unique(1)
         hi = unique(nu)
         do r = 1, ntry
            cut = lo + rng%uniform() * (hi - lo)
            call survival_split_score(values, time, status, indices, weights, cut, .false., options, min_bucket, &
               score_l, p_l, ok_l)
            call survival_split_score(values, time, status, indices, weights, cut, .true., options, min_bucket, &
               score_r, p_r, ok_r)
            call choose_numeric(score_l, p_l, ok_l, score_r, p_r, ok_r, cut, best_score, best_cut, best_nan_right, &
               best_pvalue, found)
         end do
      else
         do i = 1, nu - 1
            cut = midpoint_safe(unique(i), unique(i + 1))
            call survival_split_score(values, time, status, indices, weights, cut, .false., options, min_bucket, &
               score_l, p_l, ok_l)
            call survival_split_score(values, time, status, indices, weights, cut, .true., options, min_bucket, &
               score_r, p_r, ok_r)
            call choose_numeric(score_l, p_l, ok_l, score_r, p_r, ok_r, cut, best_score, best_cut, best_nan_right, &
               best_pvalue, found)
         end do
      end if
   end subroutine best_numeric_survival

   subroutine best_categorical_survival(values, time, status, indices, weights, ncat, options, min_bucket, rng, &
      best_score, best_mask, best_nan_right, best_pvalue, found)
      real(dp), intent(in) :: values(:), time(:), weights(:)
      integer, intent(in) :: status(:), indices(:), ncat, min_bucket
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      real(dp), intent(out) :: best_score, best_pvalue
      logical, intent(out) :: best_mask(:), best_nan_right, found
      logical, allocatable :: mask(:)
      integer, allocatable :: present(:)
      integer :: r, ntry, npresent
      integer(i64) :: split_id, max_split_id
      real(dp) :: score_l, score_r, p_l, p_r
      logical :: ok_l, ok_r

      allocate(mask(ncat))
      best_score = -huge(1.0_dp)
      best_mask = .false.
      best_nan_right = .false.
      best_pvalue = 1.0_dp
      found = .false.

      if (options%split_rule == RANGER_SPLIT_EXTRATREES) then
         ntry = max(1, options%num_random_splits)
         do r = 1, ntry
            call random_extratrees_factor_mask(rng, values, indices, ncat, mask)
            call categorical_survival_score(values, time, status, indices, weights, mask, .false., options, &
               min_bucket, score_l, p_l, ok_l)
            call categorical_survival_score(values, time, status, indices, weights, mask, .true., options, &
               min_bucket, score_r, p_r, ok_r)
            call choose_mask(score_l, p_l, ok_l, score_r, p_r, ok_r, mask, best_score, best_mask, &
               best_nan_right, best_pvalue, found)
         end do
      else
         call collect_present_categories(values, indices, ncat, present, npresent)
         if (npresent < 2) return
         max_split_id = shiftl(1_i64, npresent - 1) - 1_i64
         do split_id = 1_i64, max_split_id
            call partition_mask_from_id(present, npresent, split_id, mask)
            call categorical_survival_score(values, time, status, indices, weights, mask, .false., options, &
               min_bucket, score_l, p_l, ok_l)
            call categorical_survival_score(values, time, status, indices, weights, mask, .true., options, &
               min_bucket, score_r, p_r, ok_r)
            call choose_mask(score_l, p_l, ok_l, score_r, p_r, ok_r, mask, best_score, best_mask, &
               best_nan_right, best_pvalue, found)
         end do
      end if
   end subroutine best_categorical_survival

   subroutine survival_split_score(values, time, status, indices, weights, cut, nan_right, options, min_bucket, &
      score, pvalue, ok)
      real(dp), intent(in) :: values(:), time(:), weights(:), cut
      integer, intent(in) :: status(:), indices(:), min_bucket
      logical, intent(in) :: nan_right
      type(ranger_options), intent(in) :: options
      real(dp), intent(out) :: score, pvalue
      logical, intent(out) :: ok
      logical, allocatable :: left(:)
      integer :: i
      allocate(left(size(indices)))
      do i = 1, size(indices)
         if (ieee_is_nan(values(indices(i)))) then
            left(i) = .not. nan_right
         else
            left(i) = values(indices(i)) <= cut
         end if
      end do
      call survival_partition_score(values, time, status, indices, weights, left, options, min_bucket, score, pvalue, ok)
   end subroutine survival_split_score

   subroutine categorical_survival_score(values, time, status, indices, weights, mask, nan_right, options, &
      min_bucket, score, pvalue, ok)
      real(dp), intent(in) :: values(:), time(:), weights(:)
      integer, intent(in) :: status(:), indices(:), min_bucket
      logical, intent(in) :: mask(:), nan_right
      type(ranger_options), intent(in) :: options
      real(dp), intent(out) :: score, pvalue
      logical, intent(out) :: ok
      logical, allocatable :: left(:)
      integer :: i, c
      allocate(left(size(indices)))
      do i = 1, size(indices)
         if (ieee_is_nan(values(indices(i)))) then
            left(i) = .not. nan_right
         else
            c = nint(values(indices(i)))
            if (c >= 1 .and. c <= size(mask)) then
               left(i) = mask(c)
            else
               left(i) = .false.
            end if
         end if
      end do
      call survival_partition_score(values, time, status, indices, weights, left, options, min_bucket, score, pvalue, ok)
   end subroutine categorical_survival_score

   subroutine survival_partition_score(values, time, status, indices, weights, left, options, min_bucket, &
      score, pvalue, ok)
      real(dp), intent(in) :: values(:), time(:), weights(:)
      integer, intent(in) :: status(:), indices(:), min_bucket
      logical, intent(in) :: left(:)
      type(ranger_options), intent(in) :: options
      real(dp), intent(out) :: score, pvalue
      logical, intent(out) :: ok
      real(dp) :: wl, wr, prop

      wl = sum(pack(weights, left))
      wr = sum(pack(weights, .not. left))
      ok = wl >= real(min_bucket, dp) .and. wr >= real(min_bucket, dp)
      if (.not. ok) then
         score = -huge(1.0_dp)
         pvalue = 1.0_dp
         return
      end if
      select case (options%split_rule)
      case (RANGER_SPLIT_AUC, RANGER_SPLIT_AUC_IGNORE_TIES)
         score = auc_split_score(values, time, status, indices, left, options%split_rule == RANGER_SPLIT_AUC_IGNORE_TIES)
         pvalue = 0.0_dp
      case (RANGER_SPLIT_MAXSTAT)
         prop = wl / (wl + wr)
         if (prop < options%minprop .or. prop > 1.0_dp - options%minprop) then
            ok = .false.
            score = -huge(1.0_dp)
            pvalue = 1.0_dp
            return
         end if
         score = logrank_split_score(time, status, indices, weights, left)
         pvalue = erfc(score / sqrt(2.0_dp))
      case default
         score = logrank_split_score(time, status, indices, weights, left)
         pvalue = 0.0_dp
      end select
   end subroutine survival_partition_score

   real(dp) function logrank_split_score(time, status, indices, weights, left) result(score)
      real(dp), intent(in) :: time(:), weights(:)
      integer, intent(in) :: status(:), indices(:)
      logical, intent(in) :: left(:)
      real(dp), allocatable :: event_times(:)
      real(dp) :: numerator, denom2, yi, yi1, di, di1
      integer :: i, j, ne, nt

      allocate(event_times(size(indices)))
      ne = 0
      do i = 1, size(indices)
         if (status(indices(i)) /= 1) cycle
         if (ne == 0) then
            ne = 1
            event_times(ne) = time(indices(i))
         else if (.not. any(abs(event_times(1:ne) - time(indices(i))) <= 0.0_dp)) then
            ne = ne + 1
            event_times(ne) = time(indices(i))
         end if
      end do
      call sort_real(event_times, ne)
      numerator = 0.0_dp
      denom2 = 0.0_dp
      do nt = 1, ne
         yi = 0.0_dp
         yi1 = 0.0_dp
         di = 0.0_dp
         di1 = 0.0_dp
         do j = 1, size(indices)
            if (time(indices(j)) >= event_times(nt)) then
               yi = yi + weights(j)
               if (.not. left(j)) yi1 = yi1 + weights(j)
            end if
            if (status(indices(j)) == 1 .and. abs(time(indices(j)) - event_times(nt)) <= 0.0_dp) then
               di = di + weights(j)
               if (.not. left(j)) di1 = di1 + weights(j)
            end if
         end do
         if (yi < 2.0_dp .or. yi1 <= 0.0_dp) cycle
         if (di > 0.0_dp) then
            numerator = numerator + di1 - yi1 * di / yi
            denom2 = denom2 + (yi1 / yi) * (1.0_dp - yi1 / yi) * ((yi - di) / (yi - 1.0_dp)) * di
         end if
      end do
      if (denom2 > 0.0_dp) then
         score = abs(numerator / sqrt(denom2))
      else
         score = -huge(1.0_dp)
      end if
   end function logrank_split_score

   real(dp) function auc_split_score(values, time, status, indices, left, ignore_ties) result(score)
      real(dp), intent(in) :: values(:), time(:)
      integer, intent(in) :: status(:), indices(:)
      logical, intent(in) :: left(:), ignore_ties
      integer :: i, j
      real(dp) :: balance, total
      logical :: comparable, do_nothing, earlier_left

      balance = 0.0_dp
      total = 0.0_dp
      do i = 1, size(indices) - 1
         do j = i + 1, size(indices)
            comparable = .true.
            do_nothing = .false.
            if (time(indices(i)) < time(indices(j))) then
               if (status(indices(i)) == 0) comparable = .false.
               earlier_left = left(i)
            else if (time(indices(j)) < time(indices(i))) then
               if (status(indices(j)) == 0) comparable = .false.
               earlier_left = left(j)
            else
               if (status(indices(i)) == 0 .or. status(indices(j)) == 0) then
                  comparable = .false.
               else if (ignore_ties) then
                  comparable = .false.
               else if (abs(values(indices(i)) - values(indices(j))) <= 0.0_dp) then
                  comparable = .false.
               else
                  do_nothing = .true.
                  earlier_left = .false.
               end if
            end if
            if (.not. comparable) cycle
            total = total + 1.0_dp
            if (do_nothing .or. left(i) .eqv. left(j)) cycle
            if (earlier_left) then
               balance = balance + 1.0_dp
            else
               balance = balance - 1.0_dp
            end if
         end do
      end do
      if (total > 0.0_dp) then
         score = abs(balance) / (2.0_dp * total)
      else
         score = -huge(1.0_dp)
      end if
   end function auc_split_score

   subroutine sort_real(x, n)
      real(dp), intent(inout) :: x(:)
      integer, intent(in) :: n
      integer :: i, j
      real(dp) :: key
      do i = 2, n
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine sort_real

   subroutine choose_numeric(score_l, p_l, ok_l, score_r, p_r, ok_r, cut, best_score, best_cut, best_nan_right, &
      best_p, found)
      real(dp), intent(in) :: score_l, p_l, score_r, p_r, cut
      logical, intent(in) :: ok_l, ok_r
      real(dp), intent(inout) :: best_score, best_cut, best_p
      logical, intent(inout) :: best_nan_right, found
      if (ok_l .and. (.not. found .or. score_l > best_score)) then
         found = .true.
         best_score = score_l
         best_cut = cut
         best_nan_right = .false.
         best_p = p_l
      end if
      if (ok_r .and. (.not. found .or. score_r > best_score)) then
         found = .true.
         best_score = score_r
         best_cut = cut
         best_nan_right = .true.
         best_p = p_r
      end if
   end subroutine choose_numeric

   subroutine choose_mask(score_l, p_l, ok_l, score_r, p_r, ok_r, mask, best_score, best_mask, best_nan_right, &
      best_p, found)
      real(dp), intent(in) :: score_l, p_l, score_r, p_r
      logical, intent(in) :: ok_l, ok_r, mask(:)
      real(dp), intent(inout) :: best_score, best_p
      logical, intent(inout) :: best_mask(:), best_nan_right, found
      if (ok_l .and. (.not. found .or. score_l > best_score)) then
         found = .true.
         best_score = score_l
         best_mask = .false.
         best_mask(1:size(mask)) = mask
         best_nan_right = .false.
         best_p = p_l
      end if
      if (ok_r .and. (.not. found .or. score_r > best_score)) then
         found = .true.
         best_score = score_r
         best_mask = .false.
         best_mask(1:size(mask)) = mask
         best_nan_right = .true.
         best_p = p_r
      end if
   end subroutine choose_mask

   subroutine split_indices(values, indices, weights, ncat, cut, mask, nan_right, left_idx, left_w, nleft, &
      right_idx, right_w, nright)
      real(dp), intent(in) :: values(:), weights(:), cut
      integer, intent(in) :: indices(:), ncat
      logical, intent(in) :: mask(:), nan_right
      integer, intent(out) :: left_idx(:), right_idx(:), nleft, nright
      real(dp), intent(out) :: left_w(:), right_w(:)
      integer :: i
      nleft = 0
      nright = 0
      do i = 1, size(indices)
         if (go_left_value(values(indices(i)), ncat, cut, mask, nan_right)) then
            nleft = nleft + 1
            left_idx(nleft) = indices(i)
            left_w(nleft) = weights(i)
         else
            nright = nright + 1
            right_idx(nright) = indices(i)
            right_w(nright) = weights(i)
         end if
      end do
   end subroutine split_indices

   subroutine predict_survival_tree(tree, x, ncat, chf, terminal_node)
      type(ranger_tree), intent(in) :: tree
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: ncat(:)
      real(dp), intent(out), optional :: chf(:,:)
      integer, intent(out), optional :: terminal_node(:)
      integer :: i, node, v
      do i = 1, size(x, 1)
         node = 1
         do while (tree%status(node) == RANGER_INTERIOR)
            v = tree%split_var(node)
            if (go_left_value(x(i, v), ncat(v), tree%split_value(node), tree%cat_left(:, node), &
               tree%nan_go_right(node))) then
               node = tree%left(node)
            else
               node = tree%right(node)
            end if
         end do
         if (present(chf)) chf(i, :) = tree%chf(:, node)
         if (present(terminal_node)) terminal_node(i) = node
      end do
   end subroutine predict_survival_tree

   real(dp) function concordance_index(time, status, risk) result(cindex)
      real(dp), intent(in) :: time(:), risk(:)
      integer, intent(in) :: status(:)
      integer :: i, j
      real(dp) :: concordant, comparable, contribution

      concordant = 0.0_dp
      comparable = 0.0_dp
      do i = 1, size(time) - 1
         do j = i + 1, size(time)
            if (.not. permissible_pair(time(i), status(i), time(j), status(j))) cycle
            contribution = pair_concordance(time(i), time(j), risk(i), risk(j))
            concordant = concordant + contribution
            comparable = comparable + 1.0_dp
         end do
      end do
      if (comparable > 0.0_dp) then
         cindex = concordant / comparable
      else
         cindex = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
   end function concordance_index

   subroutine concordance_casewise(time, status, risk, prediction_error)
      real(dp), intent(in) :: time(:), risk(:)
      integer, intent(in) :: status(:)
      real(dp), intent(out) :: prediction_error(:)
      real(dp), allocatable :: concordant(:), permissible(:)
      real(dp) :: contribution
      integer :: i, j

      if (size(status) /= size(time) .or. size(risk) /= size(time) .or. size(prediction_error) /= size(time)) &
         error stop 'concordance_casewise: incompatible dimensions'
      allocate(concordant(size(time)), permissible(size(time)))
      concordant = 0.0_dp
      permissible = 0.0_dp
      do i = 1, size(time) - 1
         do j = i + 1, size(time)
            if (.not. permissible_pair(time(i), status(i), time(j), status(j))) cycle
            contribution = pair_concordance(time(i), time(j), risk(i), risk(j))
            concordant(i) = concordant(i) + contribution
            concordant(j) = concordant(j) + contribution
            permissible(i) = permissible(i) + 1.0_dp
            permissible(j) = permissible(j) + 1.0_dp
         end do
      end do
      do i = 1, size(time)
         if (permissible(i) > 0.0_dp) then
            prediction_error(i) = 1.0_dp - concordant(i) / permissible(i)
         else
            prediction_error(i) = ieee_value(0.0_dp, ieee_quiet_nan)
         end if
      end do
   end subroutine concordance_casewise

   pure logical function permissible_pair(time_i, status_i, time_j, status_j) result(permissible)
      real(dp), intent(in) :: time_i, time_j
      integer, intent(in) :: status_i, status_j

      permissible = .true.
      if (time_i < time_j .and. status_i == 0) permissible = .false.
      if (time_j < time_i .and. status_j == 0) permissible = .false.
      if (time_i <= time_j .and. time_i >= time_j .and. status_i == status_j) permissible = .false.
   end function permissible_pair

   pure real(dp) function pair_concordance(time_i, time_j, risk_i, risk_j) result(value)
      real(dp), intent(in) :: time_i, time_j, risk_i, risk_j

      value = 0.0_dp
      if (time_i < time_j .and. risk_i > risk_j) then
         value = 1.0_dp
      else if (time_j < time_i .and. risk_j > risk_i) then
         value = 1.0_dp
      else if (abs(risk_i - risk_j) <= 0.0_dp) then
         value = 0.5_dp
      end if
   end function pair_concordance

end module ranger_tree_survival
