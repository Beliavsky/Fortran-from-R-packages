! SPDX-License-Identifier: GPL-3.0-only
! Modern Fortran translation of ranger 0.18.0 computational code.
! Upstream C++ core: Copyright (c) 2014-2018 Marvin N. Wright; MIT licensed.
! The upstream R package is GPL-3. See NOTICE.md for full provenance.
module ranger_tree_regression
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
   use r_kinds, only : dp, i64
   use ranger_rng, only : ranger_rng_state
   use ranger_types, only : ranger_options, ranger_tree, RANGER_INTERIOR
   use ranger_types, only : RANGER_SPLIT_STANDARD, RANGER_SPLIT_MAXSTAT, RANGER_SPLIT_EXTRATREES
   use ranger_types, only : RANGER_SPLIT_BETA, RANGER_SPLIT_POISSON
   use ranger_types, only : RANGER_UNORDERED_PARTITION
   use ranger_tree_common, only : initialize_tree, go_left_value, candidate_variables, collect_sorted_unique
   use ranger_tree_common, only : midpoint_safe, regularized_score
   use ranger_tree_common, only : collect_present_categories, partition_mask_from_id
   use ranger_tree_common, only : random_extratrees_factor_mask
   use ranger_tree_common, only : base_predictor_count, unpermuted_var_id, is_shadow_variable
   use ranger_maxstat, only : average_ranks, maxstat_best, maxstat_pvalue_lau92, maxstat_pvalue_lau94
   use ranger_maxstat, only : adjust_pvalues_bh
   implicit none
   private

   public :: build_regression_tree, predict_regression_tree

contains

   subroutine build_regression_tree(x, y, obs_index, obs_weight, ncat, options, rng, used_global, tree)
      real(dp), intent(in) :: x(:,:), y(:), obs_weight(:)
      integer, intent(in) :: obs_index(:), ncat(:)
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      logical, intent(inout) :: used_global(:)
      type(ranger_tree), intent(out) :: tree
      integer :: maxnodes, maxcat, next_node, min_node, min_bucket, pbase
      real(dp) :: root_mean

      pbase = base_predictor_count(size(x, 2), options)
      maxnodes = max(3, 2 * size(obs_index) + 1)
      maxcat = max(1, maxval(ncat))
      min_node = options%min_node_size
      if (min_node <= 0) min_node = 5
      min_bucket = options%min_bucket
      if (min_bucket <= 0) min_bucket = 1
      call initialize_tree(tree, maxnodes, maxcat, pbase)
      next_node = 1
      root_mean = weighted_mean(y, obs_index, obs_weight)
      call grow_node(obs_index, obs_weight, 1, 0, root_mean)
      tree%n_nodes = next_node

   contains

      recursive subroutine grow_node(indices, weights, node, depth, parent_mean)
         integer, intent(in) :: indices(:), node, depth
         real(dp), intent(in) :: weights(:), parent_mean
         integer, allocatable :: left_idx(:), right_idx(:)
         real(dp), allocatable :: left_w(:), right_w(:)
         logical, allocatable :: best_mask(:)
         integer :: best_var, base_var, nleft, nright, left_node, right_node
         real(dp) :: best_cut, best_score, best_importance, node_mean, node_sum, weight_sum
         real(dp) :: alpha_smooth, importance_sign
         logical :: best_nan_right, found

         weight_sum = sum(weights)
         node_sum = weighted_sum(y, indices, weights)
         if (weight_sum > 0.0_dp) then
            node_mean = node_sum / weight_sum
         else
            node_mean = parent_mean
         end if
         if (options%split_rule == RANGER_SPLIT_POISSON .and. node_sum <= 0.0_dp .and. node > 1) then
            alpha_smooth = weight_sum * max(parent_mean, 0.0_dp) / &
               (weight_sum * max(parent_mean, 0.0_dp) + options%poisson_tau)
            node_mean = alpha_smooth * node_mean + (1.0_dp - alpha_smooth) * parent_mean
         end if
         tree%node_n(node) = nint(weight_sum)
         tree%node_mean(node) = node_mean

         if (size(indices) <= min_node) return
         if (options%max_depth > 0 .and. depth >= options%max_depth) return
         if (maxval(y(indices)) - minval(y(indices)) <= 0.0_dp) return
         if (next_node + 2 > size(tree%left)) return

         allocate(best_mask(maxcat))
         call best_reg_split(x, y, indices, weights, ncat, options, min_bucket, depth, rng, used_global, best_var, &
            best_cut, best_mask, best_nan_right, best_score, best_importance, found)
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
         tree%impurity_decrease(base_var) = tree%impurity_decrease(base_var) + importance_sign * best_importance
         used_global(base_var) = .true.

         left_node = next_node + 1
         right_node = next_node + 2
         next_node = next_node + 2
         tree%left(node) = left_node
         tree%right(node) = right_node
         call grow_node(left_idx(1:nleft), left_w(1:nleft), left_node, depth + 1, node_mean)
         call grow_node(right_idx(1:nright), right_w(1:nright), right_node, depth + 1, node_mean)
      end subroutine grow_node

   end subroutine build_regression_tree

   subroutine best_reg_split(x, y, indices, weights, ncat, options, min_bucket, depth, rng, used_global, &
      best_var, best_cut, best_mask, best_nan_right, best_stat, best_importance, found)
      real(dp), intent(in) :: x(:,:), y(:), weights(:)
      integer, intent(in) :: indices(:), ncat(:), min_bucket, depth
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      logical, intent(in) :: used_global(:)
      integer, intent(out) :: best_var
      real(dp), intent(out) :: best_cut, best_stat, best_importance
      logical, intent(out) :: best_mask(:), best_nan_right, found
      integer, allocatable :: vars(:)
      logical, allocatable :: mask(:)
      integer :: i, v, base_v, nvars, mtry, pbase
      real(dp) :: score, cut, parent_impurity, adjusted, parent_adjusted, candidate_importance
      real(dp) :: pvalue, best_p
      logical :: nan_right, ok, negative

      if (options%split_rule == RANGER_SPLIT_MAXSTAT) then
         call best_reg_split_maxstat(x, y, indices, ncat, options, rng, best_var, best_cut, best_mask, &
            best_nan_right, best_stat, found)
         best_importance = best_stat
         return
      end if

      parent_impurity = weighted_sum(y, indices, weights) ** 2 / sum(weights)
      pbase = base_predictor_count(size(x, 2), options)
      mtry = options%mtry
      if (mtry <= 0) mtry = max(1, int(sqrt(real(pbase, dp))))
      call candidate_variables(rng, size(x, 2), mtry, options, vars, nvars)
      allocate(mask(size(best_mask)))

      best_var = 0
      best_cut = 0.0_dp
      best_stat = -huge(1.0_dp)
      best_importance = 0.0_dp
      best_mask = .false.
      best_nan_right = .false.
      best_p = 1.0_dp
      found = .false.
      do i = 1, nvars
         v = vars(i)
         mask = .false.
         call best_reg_variable(x(:, v), y, indices, weights, ncat(v), options, min_bucket, rng, score, cut, mask, &
            nan_right, pvalue, ok)
         if (.not. ok) cycle
         if (options%split_rule == RANGER_SPLIT_BETA) then
            negative = .true.
         else if (options%split_rule == RANGER_SPLIT_POISSON) then
            negative = score <= 0.0_dp
         else
            negative = .false.
         end if
         base_v = unpermuted_var_id(v, size(x, 2), options)
         adjusted = regularized_score(score, base_v, depth, options, used_global, negative)
         parent_adjusted = regularized_score(parent_impurity, base_v, depth, options, used_global)
         candidate_importance = adjusted - parent_adjusted
         if (.not. found .or. adjusted > best_stat) then
            found = .true.
            best_var = v
            best_cut = cut
            best_mask = mask
            best_nan_right = nan_right
            best_stat = adjusted
            best_importance = candidate_importance
            best_p = pvalue
         end if
      end do
   end subroutine best_reg_split

   subroutine best_reg_split_maxstat(x, y, indices, ncat, options, rng, best_var, best_cut, best_mask, &
      best_nan_right, best_stat, found)
      real(dp), intent(in) :: x(:,:), y(:)
      integer, intent(in) :: indices(:), ncat(:)
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      integer, intent(out) :: best_var
      real(dp), intent(out) :: best_cut, best_stat
      logical, intent(out) :: best_mask(:), best_nan_right, found
      integer, allocatable :: vars(:), num_left(:)
      real(dp), allocatable :: response(:), ranks(:), xnode(:), pvalue(:), adjusted(:), stat(:), cut(:)
      logical, allocatable :: valid(:)
      integer :: nvars, mtry, pbase, i, v, nvalid, imin, c
      real(dp) :: p92, p94
      logical :: ok

      pbase = base_predictor_count(size(x, 2), options)
      mtry = options%mtry
      if (mtry <= 0) mtry = max(1, int(sqrt(real(pbase, dp))))
      call candidate_variables(rng, size(x, 2), mtry, options, vars, nvars)
      allocate(response(size(indices)), ranks(size(indices)), xnode(size(indices)))
      allocate(pvalue(nvars), adjusted(nvars), stat(nvars), cut(nvars), valid(nvars))
      response = y(indices)
      call average_ranks(response, ranks)
      pvalue = 1.0_dp
      adjusted = 1.0_dp
      stat = -1.0_dp
      cut = 0.0_dp
      valid = .false.

      do i = 1, nvars
         v = vars(i)
         xnode = x(indices, v)
         call maxstat_best(ranks, xnode, options%minprop, 1.0_dp - options%minprop, stat(i), cut(i), num_left, ok)
         if (.not. ok) cycle
         p92 = maxstat_pvalue_lau92(stat(i), options%minprop, 1.0_dp - options%minprop)
         p94 = maxstat_pvalue_lau94(stat(i), options%minprop, 1.0_dp - options%minprop, size(indices), num_left)
         pvalue(i) = min(p92, p94)
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
      call adjust_valid_pvalues(pvalue, valid, adjusted)

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
   end subroutine best_reg_split_maxstat

   subroutine adjust_valid_pvalues(pvalue, valid, adjusted)
      real(dp), intent(in) :: pvalue(:)
      logical, intent(in) :: valid(:)
      real(dp), intent(out) :: adjusted(size(pvalue))
      real(dp), allocatable :: pv(:), av(:)
      integer, allocatable :: map(:)
      integer :: i, j, n

      n = count(valid)
      allocate(pv(n), av(n), map(n))
      j = 0
      do i = 1, size(pvalue)
         if (.not. valid(i)) cycle
         j = j + 1
         pv(j) = pvalue(i)
         map(j) = i
      end do
      call adjust_pvalues_bh(pv, av)
      adjusted = 1.0_dp
      do j = 1, n
         adjusted(map(j)) = av(j)
      end do
   end subroutine adjust_valid_pvalues

   subroutine best_reg_variable(values, y, indices, weights, ncat, options, min_bucket, rng, best_score, best_cut, &
      best_mask, best_nan_right, best_pvalue, found)
      real(dp), intent(in) :: values(:), y(:), weights(:)
      integer, intent(in) :: indices(:), ncat, min_bucket
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      real(dp), intent(out) :: best_score, best_cut, best_pvalue
      logical, intent(out) :: best_mask(:), best_nan_right, found

      if (ncat <= 1 .or. options%respect_unordered_factors /= RANGER_UNORDERED_PARTITION) then
         call best_numeric_reg(values, y, indices, weights, options, min_bucket, rng, best_score, best_cut, &
            best_nan_right, best_pvalue, found)
         best_mask = .false.
      else
         call best_categorical_reg(values, y, indices, weights, ncat, options, min_bucket, rng, best_score, &
            best_mask, best_nan_right, best_pvalue, found)
         best_cut = 0.0_dp
      end if
   end subroutine best_reg_variable

   subroutine best_numeric_reg(values, y, indices, weights, options, min_bucket, rng, best_score, best_cut, &
      best_nan_right, best_pvalue, found)
      real(dp), intent(in) :: values(:), y(:), weights(:)
      integer, intent(in) :: indices(:), min_bucket
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
            call reg_split_score(values, y, indices, weights, cut, .false., options, min_bucket, score_l, p_l, ok_l)
            call reg_split_score(values, y, indices, weights, cut, .true., options, min_bucket, score_r, p_r, ok_r)
            call choose_numeric(score_l, p_l, ok_l, score_r, p_r, ok_r, cut, best_score, best_cut, best_nan_right, &
               best_pvalue, found)
         end do
      else
         do i = 1, nu - 1
            cut = midpoint_safe(unique(i), unique(i + 1))
            call reg_split_score(values, y, indices, weights, cut, .false., options, min_bucket, score_l, p_l, ok_l)
            call reg_split_score(values, y, indices, weights, cut, .true., options, min_bucket, score_r, p_r, ok_r)
            call choose_numeric(score_l, p_l, ok_l, score_r, p_r, ok_r, cut, best_score, best_cut, best_nan_right, &
               best_pvalue, found)
         end do
      end if
   end subroutine best_numeric_reg

   subroutine best_categorical_reg(values, y, indices, weights, ncat, options, min_bucket, rng, best_score, &
      best_mask, best_nan_right, best_pvalue, found)
      real(dp), intent(in) :: values(:), y(:), weights(:)
      integer, intent(in) :: indices(:), ncat, min_bucket
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
            call categorical_reg_score(values, y, indices, weights, mask, .false., options, min_bucket, &
               score_l, p_l, ok_l)
            call categorical_reg_score(values, y, indices, weights, mask, .true., options, min_bucket, &
               score_r, p_r, ok_r)
            call choose_mask(score_l, p_l, ok_l, score_r, p_r, ok_r, mask, best_score, best_mask, &
               best_nan_right, best_pvalue, found)
         end do
      else
         call collect_present_categories(values, indices, ncat, present, npresent)
         if (npresent < 2) return
         max_split_id = shiftl(1_i64, npresent - 1) - 1_i64
         do split_id = 1_i64, max_split_id
            call partition_mask_from_id(present, npresent, split_id, mask)
            call categorical_reg_score(values, y, indices, weights, mask, .false., options, min_bucket, &
               score_l, p_l, ok_l)
            call categorical_reg_score(values, y, indices, weights, mask, .true., options, min_bucket, &
               score_r, p_r, ok_r)
            call choose_mask(score_l, p_l, ok_l, score_r, p_r, ok_r, mask, best_score, best_mask, &
               best_nan_right, best_pvalue, found)
         end do
      end if
   end subroutine best_categorical_reg

   subroutine reg_split_score(values, y, indices, weights, cut, nan_right, options, min_bucket, score, pvalue, ok)
      real(dp), intent(in) :: values(:), y(:), weights(:), cut
      integer, intent(in) :: indices(:), min_bucket
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
      call score_partition(y, indices, weights, left, options, min_bucket, score, pvalue, ok)
   end subroutine reg_split_score

   subroutine categorical_reg_score(values, y, indices, weights, mask, nan_right, options, min_bucket, &
      score, pvalue, ok)
      real(dp), intent(in) :: values(:), y(:), weights(:)
      integer, intent(in) :: indices(:), min_bucket
      logical, intent(in) :: mask(:), nan_right
      type(ranger_options), intent(in) :: options
      real(dp), intent(out) :: score, pvalue
      logical, intent(out) :: ok
      logical, allocatable :: left(:)
      integer :: i, cat

      allocate(left(size(indices)))
      do i = 1, size(indices)
         if (ieee_is_nan(values(indices(i)))) then
            left(i) = .not. nan_right
         else
            cat = nint(values(indices(i)))
            if (cat >= 1 .and. cat <= size(mask)) then
               left(i) = mask(cat)
            else
               left(i) = .false.
            end if
         end if
      end do
      call score_partition(y, indices, weights, left, options, min_bucket, score, pvalue, ok)
   end subroutine categorical_reg_score

   subroutine score_partition(y, indices, weights, left, options, min_bucket, score, pvalue, ok)
      real(dp), intent(in) :: y(:), weights(:)
      integer, intent(in) :: indices(:), min_bucket
      logical, intent(in) :: left(:)
      type(ranger_options), intent(in) :: options
      real(dp), intent(out) :: score, pvalue
      logical, intent(out) :: ok
      real(dp) :: wl, wr, sl, sr, ml, mr, vl, vr, phi_l, phi_r, total_w, prop
      integer :: i, nl, nr

      wl = 0.0_dp
      wr = 0.0_dp
      sl = 0.0_dp
      sr = 0.0_dp
      nl = 0
      nr = 0
      do i = 1, size(indices)
         if (left(i)) then
            wl = wl + weights(i)
            sl = sl + weights(i) * y(indices(i))
            nl = nl + nint(weights(i))
         else
            wr = wr + weights(i)
            sr = sr + weights(i) * y(indices(i))
            nr = nr + nint(weights(i))
         end if
      end do
      ok = wl >= real(min_bucket, dp) .and. wr >= real(min_bucket, dp) .and. nl > 0 .and. nr > 0
      if (.not. ok) then
         score = -huge(1.0_dp)
         pvalue = 1.0_dp
         return
      end if
      ml = sl / wl
      mr = sr / wr
      pvalue = 0.0_dp

      select case (options%split_rule)
      case (RANGER_SPLIT_BETA)
         if (nl < 2 .or. nr < 2) then
            ok = .false.
            score = -huge(1.0_dp)
            return
         end if
         vl = weighted_variance_partition(y, indices, weights, left, .true., ml)
         vr = weighted_variance_partition(y, indices, weights, left, .false., mr)
         if (vl < epsilon(1.0_dp) .or. vr < epsilon(1.0_dp)) then
            ok = .false.
            score = -huge(1.0_dp)
            return
         end if
         phi_l = ml * (1.0_dp - ml) / vl - 1.0_dp
         phi_r = mr * (1.0_dp - mr) / vr - 1.0_dp
         score = beta_partition_loglik(y, indices, weights, left, ml, mr, phi_l, phi_r)
      case (RANGER_SPLIT_POISSON)
         score = xlogy(sl, ml) + xlogy(sr, mr)
      case (RANGER_SPLIT_MAXSTAT)
         total_w = wl + wr
         prop = wl / total_w
         if (prop < options%minprop .or. prop > 1.0_dp - options%minprop) then
            ok = .false.
            score = -huge(1.0_dp)
            pvalue = 1.0_dp
            return
         end if
         vl = weighted_variance_all(y, indices, weights)
         if (vl <= 0.0_dp) then
            ok = .false.
            score = -huge(1.0_dp)
            pvalue = 1.0_dp
            return
         end if
         score = abs(ml - mr) / sqrt(vl * (1.0_dp / wl + 1.0_dp / wr))
         pvalue = erfc(score / sqrt(2.0_dp))
      case default
         score = sl * sl / wl + sr * sr / wr
      end select
      if (.not. ieee_is_finite(score)) ok = .false.
   end subroutine score_partition


   pure real(dp) function xlogy(x, yval) result(value)
      real(dp), intent(in) :: x, yval
      if (x <= 0.0_dp) then
         value = 0.0_dp
      else if (yval <= 0.0_dp) then
         value = -huge(1.0_dp)
      else
         value = x * log(yval)
      end if
   end function xlogy

   pure real(dp) function beta_loglik_one(yval, mean, phi) result(value)
      real(dp), intent(in) :: yval, mean, phi
      real(dp) :: yy, mm, pp, eps

      eps = epsilon(1.0_dp)
      yy = yval
      if (yy < eps) then
         yy = eps
      else if (yy >= 1.0_dp) then
         yy = 1.0_dp - eps
      end if
      mm = mean
      if (mm < eps) then
         mm = eps
      else if (mm >= 1.0_dp) then
         mm = 1.0_dp - eps
      end if
      pp = phi
      if (pp < eps) pp = eps
      value = log_gamma(pp) - log_gamma(mm * pp) - log_gamma((1.0_dp - mm) * pp) + &
         (mm * pp - 1.0_dp) * log(yy) + ((1.0_dp - mm) * pp - 1.0_dp) * log(1.0_dp - yy)
   end function beta_loglik_one

   real(dp) function beta_partition_loglik(y, indices, weights, left, ml, mr, phi_l, phi_r) result(value)
      real(dp), intent(in) :: y(:), weights(:), ml, mr, phi_l, phi_r
      integer, intent(in) :: indices(:)
      logical, intent(in) :: left(:)
      integer :: i
      value = 0.0_dp
      do i = 1, size(indices)
         if (left(i)) then
            value = value + weights(i) * beta_loglik_one(y(indices(i)), ml, phi_l)
         else
            value = value + weights(i) * beta_loglik_one(y(indices(i)), mr, phi_r)
         end if
      end do
   end function beta_partition_loglik


   real(dp) function weighted_variance_partition(y, indices, weights, left, select_left, mean) result(v)
      real(dp), intent(in) :: y(:), weights(:), mean
      integer, intent(in) :: indices(:)
      logical, intent(in) :: left(:), select_left
      real(dp) :: w
      integer :: i
      v = 0.0_dp
      w = 0.0_dp
      do i = 1, size(indices)
         if (left(i) .neqv. select_left) cycle
         v = v + weights(i) * (y(indices(i)) - mean) ** 2
         w = w + weights(i)
      end do
      if (w > 1.0_dp) then
         v = v / (w - 1.0_dp)
      else
         v = 0.0_dp
      end if
   end function weighted_variance_partition

   real(dp) function weighted_variance_all(y, indices, weights) result(v)
      real(dp), intent(in) :: y(:), weights(:)
      integer, intent(in) :: indices(:)
      real(dp) :: m, w
      integer :: i
      m = weighted_mean(y, indices, weights)
      w = sum(weights)
      v = 0.0_dp
      do i = 1, size(indices)
         v = v + weights(i) * (y(indices(i)) - m) ** 2
      end do
      if (w > 1.0_dp) then
         v = v / (w - 1.0_dp)
      else
         v = 0.0_dp
      end if
   end function weighted_variance_all

   pure real(dp) function weighted_sum(y, indices, weights) result(s)
      real(dp), intent(in) :: y(:), weights(:)
      integer, intent(in) :: indices(:)
      integer :: i
      s = 0.0_dp
      do i = 1, size(indices)
         s = s + weights(i) * y(indices(i))
      end do
   end function weighted_sum

   pure real(dp) function weighted_mean(y, indices, weights) result(m)
      real(dp), intent(in) :: y(:), weights(:)
      integer, intent(in) :: indices(:)
      real(dp) :: w
      w = sum(weights)
      if (w > 0.0_dp) then
         m = weighted_sum(y, indices, weights) / w
      else
         m = 0.0_dp
      end if
   end function weighted_mean


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

   subroutine predict_regression_tree(tree, x, ncat, prediction, terminal_node)
      type(ranger_tree), intent(in) :: tree
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: ncat(:)
      real(dp), intent(out), optional :: prediction(:)
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
         if (present(prediction)) prediction(i) = tree%node_mean(node)
         if (present(terminal_node)) terminal_node(i) = node
      end do
   end subroutine predict_regression_tree

end module ranger_tree_regression
