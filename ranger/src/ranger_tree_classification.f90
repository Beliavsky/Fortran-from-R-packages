! SPDX-License-Identifier: GPL-3.0-only
! Modern Fortran translation of ranger 0.18.0 computational code.
! Upstream C++ core: Copyright (c) 2014-2018 Marvin N. Wright; MIT licensed.
! The upstream R package is GPL-3. See NOTICE.md for full provenance.
module ranger_tree_classification
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp, i64
   use ranger_rng, only : ranger_rng_state
   use ranger_types, only : ranger_options, ranger_tree, RANGER_INTERIOR, RANGER_TERMINAL
   use ranger_types, only : RANGER_SPLIT_STANDARD, RANGER_SPLIT_EXTRATREES, RANGER_SPLIT_HELLINGER
   use ranger_types, only : RANGER_UNORDERED_PARTITION
   use ranger_tree_common, only : initialize_tree, go_left_value, candidate_variables, collect_sorted_unique
   use ranger_tree_common, only : midpoint_safe, regularized_score
   use ranger_tree_common, only : collect_present_categories, partition_mask_from_id
   use ranger_tree_common, only : random_extratrees_factor_mask
   use ranger_tree_common, only : base_predictor_count, unpermuted_var_id, is_shadow_variable
   implicit none
   private

   public :: build_classification_tree, predict_classification_tree

contains

   subroutine build_classification_tree(x, y, obs_index, obs_weight, ncat, nclass, class_weights, options, &
      min_node_class, min_bucket_class, rng, used_global, tree)
      real(dp), intent(in) :: x(:,:), obs_weight(:), class_weights(:)
      integer, intent(in) :: y(:), obs_index(:), ncat(:), nclass
      type(ranger_options), intent(in) :: options
      integer, intent(in) :: min_node_class(:), min_bucket_class(:)
      type(ranger_rng_state), intent(inout) :: rng
      logical, intent(inout) :: used_global(:)
      type(ranger_tree), intent(out) :: tree
      integer :: maxnodes, maxcat, next_node, pbase

      pbase = base_predictor_count(size(x, 2), options)
      maxnodes = max(3, 2 * size(obs_index) + 1)
      maxcat = max(1, maxval(ncat))
      call initialize_tree(tree, maxnodes, maxcat, pbase, nclass=nclass)
      next_node = 1
      call grow_node(obs_index, obs_weight, 1, 0)
      tree%n_nodes = next_node

   contains

      recursive subroutine grow_node(indices, weights, node, depth)
         integer, intent(in) :: indices(:), node, depth
         real(dp), intent(in) :: weights(:)
         real(dp), allocatable :: counts(:)
         integer, allocatable :: left_idx(:), right_idx(:)
         real(dp), allocatable :: left_w(:), right_w(:)
         logical, allocatable :: best_mask(:)
         integer :: i, best_var, base_var, nleft, nright, left_node, right_node
         real(dp) :: best_cut, best_score, best_importance, total_count, importance_sign
         logical :: best_nan_right, found

         allocate(counts(nclass))
         counts = 0.0_dp
         do i = 1, size(indices)
            counts(y(indices(i))) = counts(y(indices(i))) + weights(i)
         end do
         total_count = sum(counts)
         tree%node_n(node) = nint(total_count)
         call store_class_node(tree, node, counts, class_weights, rng)

         if (should_stop(counts, size(indices), depth, options, min_node_class)) return
         if (next_node + 2 > size(tree%left)) return

         allocate(best_mask(maxcat))
         call best_class_split(x, y, indices, weights, ncat, nclass, class_weights, options, min_bucket_class, &
            depth, rng, used_global, best_var, best_cut, best_mask, best_nan_right, best_score, &
            best_importance, found)
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
         call grow_node(left_idx(1:nleft), left_w(1:nleft), left_node, depth + 1)
         call grow_node(right_idx(1:nright), right_w(1:nright), right_node, depth + 1)
      end subroutine grow_node

   end subroutine build_classification_tree

   subroutine store_class_node(tree, node, counts, class_weights, rng)
      type(ranger_tree), intent(inout) :: tree
      integer, intent(in) :: node
      real(dp), intent(in) :: counts(:), class_weights(:)
      type(ranger_rng_state), intent(in) :: rng
      type(ranger_rng_state) :: tie_rng
      real(dp) :: total, best, value
      integer, allocatable :: tied(:)
      integer :: k, chosen, ntie

      total = sum(counts)
      if (total > 0.0_dp) tree%class_prob(:, node) = counts / total
      allocate(tied(size(counts)))
      best = 0.0_dp
      ntie = 0
      do k = 1, size(counts)
         value = class_weights(k) * counts(k)
         if (value > best) then
            best = value
            ntie = 1
            tied(1) = k
         else if (abs(value - best) <= 0.0_dp) then
            ntie = ntie + 1
            tied(ntie) = k
         end if
      end do
      if (best <= 0.0_dp .or. ntie == 0) then
         chosen = 1
      else if (ntie == 1) then
         chosen = tied(1)
      else
         ! utility.h::mostFrequentClass() takes the RNG by value, so the
         ! random tie break does not advance the tree RNG.
         tie_rng = rng
         chosen = tied(tie_rng%randint(ntie))
      end if
      tree%node_class(node) = chosen
   end subroutine store_class_node

   logical function should_stop(counts, nobs, depth, options, min_node_class) result(stop)
      real(dp), intent(in) :: counts(:)
      integer, intent(in) :: nobs, depth, min_node_class(:)
      type(ranger_options), intent(in) :: options
      integer :: k

      stop = .false.
      if (count(counts > 0.0_dp) <= 1) then
         stop = .true.
         return
      end if
      if (options%max_depth > 0 .and. depth >= options%max_depth) then
         stop = .true.
         return
      end if
      if (size(min_node_class) == 1) then
         stop = nobs <= min_node_class(1)
      else
         stop = .false.
         do k = 1, min(size(counts), size(min_node_class))
            if (counts(k) < real(min_node_class(k), dp)) then
               stop = .true.
               exit
            end if
         end do
      end if
   end function should_stop

   subroutine best_class_split(x, y, indices, weights, ncat, nclass, class_weights, options, min_bucket_class, &
      depth, rng, used_global, best_var, best_cut, best_mask, best_nan_right, best_stat, best_importance, found)
      real(dp), intent(in) :: x(:,:), weights(:), class_weights(:)
      integer, intent(in) :: y(:), indices(:), ncat(:), nclass, min_bucket_class(:), depth
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      logical, intent(in) :: used_global(:)
      integer, intent(out) :: best_var
      real(dp), intent(out) :: best_cut, best_stat, best_importance
      logical, intent(out) :: best_mask(:), best_nan_right, found
      integer, allocatable :: vars(:)
      real(dp), allocatable :: parent_counts(:)
      logical, allocatable :: mask(:)
      integer :: i, v, base_v, nvars, mtry, pbase
      real(dp) :: score, cut, parent_score, adjusted, parent_adjusted, candidate_importance
      logical :: nan_right, ok

      allocate(parent_counts(nclass), mask(size(best_mask)))
      parent_counts = 0.0_dp
      do i = 1, size(indices)
         parent_counts(y(indices(i))) = parent_counts(y(indices(i))) + weights(i)
      end do
      parent_score = gini_score(parent_counts, class_weights)
      pbase = base_predictor_count(size(x, 2), options)
      mtry = options%mtry
      if (mtry <= 0) mtry = max(1, int(sqrt(real(pbase, dp))))
      call candidate_variables(rng, size(x, 2), mtry, options, vars, nvars)

      best_var = 0
      best_cut = 0.0_dp
      best_stat = -huge(1.0_dp)
      best_importance = 0.0_dp
      best_mask = .false.
      best_nan_right = .false.
      found = .false.
      do i = 1, nvars
         v = vars(i)
         mask = .false.
         call best_class_variable(x(:, v), y, indices, weights, ncat(v), nclass, class_weights, options, &
            min_bucket_class, rng, score, cut, mask, nan_right, ok)
         if (.not. ok) cycle
         base_v = unpermuted_var_id(v, size(x, 2), options)
         adjusted = regularized_score(score, base_v, depth, options, used_global)
         if (options%split_rule == RANGER_SPLIT_HELLINGER) then
            candidate_importance = adjusted
         else
            parent_adjusted = regularized_score(parent_score, base_v, depth, options, used_global)
            candidate_importance = adjusted - parent_adjusted
         end if
         if (.not. found .or. adjusted > best_stat) then
            found = .true.
            best_var = v
            best_cut = cut
            best_mask = mask
            best_nan_right = nan_right
            best_stat = adjusted
            best_importance = candidate_importance
         end if
      end do
   end subroutine best_class_split

   subroutine best_class_variable(values, y, indices, weights, ncat, nclass, class_weights, options, &
      min_bucket_class, rng, best_score, best_cut, best_mask, best_nan_right, found)
      real(dp), intent(in) :: values(:), weights(:), class_weights(:)
      integer, intent(in) :: y(:), indices(:), ncat, nclass, min_bucket_class(:)
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      real(dp), intent(out) :: best_score, best_cut
      logical, intent(out) :: best_mask(:), best_nan_right, found

      if (ncat <= 1 .or. options%respect_unordered_factors /= RANGER_UNORDERED_PARTITION) then
         call best_numeric_class(values, y, indices, weights, nclass, class_weights, options, min_bucket_class, &
            rng, best_score, best_cut, best_nan_right, found)
         best_mask = .false.
      else
         call best_categorical_class(values, y, indices, weights, ncat, nclass, class_weights, options, &
            min_bucket_class, rng, best_score, best_mask, best_nan_right, found)
         best_cut = 0.0_dp
      end if
   end subroutine best_class_variable

   subroutine best_numeric_class(values, y, indices, weights, nclass, class_weights, options, min_bucket_class, &
      rng, best_score, best_cut, best_nan_right, found)
      real(dp), intent(in) :: values(:), weights(:), class_weights(:)
      integer, intent(in) :: y(:), indices(:), nclass, min_bucket_class(:)
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      real(dp), intent(out) :: best_score, best_cut
      logical, intent(out) :: best_nan_right, found
      real(dp), allocatable :: unique(:)
      integer :: nu, i, r, ntry
      real(dp) :: cut, lo, hi, score_l, score_r
      logical :: ok_l, ok_r

      call collect_sorted_unique(values, indices, unique, nu)
      best_score = -huge(1.0_dp)
      best_cut = 0.0_dp
      best_nan_right = .false.
      found = .false.
      if (nu < 2) return

      if (options%split_rule == RANGER_SPLIT_EXTRATREES) then
         ntry = max(1, options%num_random_splits)
         lo = unique(1)
         hi = unique(nu)
         do r = 1, ntry
            cut = lo + rng%uniform() * (hi - lo)
            call class_split_score(values, y, indices, weights, cut, .false., nclass, class_weights, options, &
               min_bucket_class, score_l, ok_l)
            call class_split_score(values, y, indices, weights, cut, .true., nclass, class_weights, options, &
               min_bucket_class, score_r, ok_r)
            call choose_nan_direction(score_l, ok_l, score_r, ok_r, cut, best_score, best_cut, best_nan_right, found)
         end do
      else
         do i = 1, nu - 1
            cut = midpoint_safe(unique(i), unique(i + 1))
            call class_split_score(values, y, indices, weights, cut, .false., nclass, class_weights, options, &
               min_bucket_class, score_l, ok_l)
            call class_split_score(values, y, indices, weights, cut, .true., nclass, class_weights, options, &
               min_bucket_class, score_r, ok_r)
            call choose_nan_direction(score_l, ok_l, score_r, ok_r, cut, best_score, best_cut, best_nan_right, found)
         end do
      end if
   end subroutine best_numeric_class

   subroutine class_split_score(values, y, indices, weights, cut, nan_right, nclass, class_weights, options, &
      min_bucket_class, score, ok)
      real(dp), intent(in) :: values(:), weights(:), cut, class_weights(:)
      integer, intent(in) :: y(:), indices(:), nclass, min_bucket_class(:)
      logical, intent(in) :: nan_right
      type(ranger_options), intent(in) :: options
      real(dp), intent(out) :: score
      logical, intent(out) :: ok
      real(dp), allocatable :: left(:), right(:)
      integer :: i

      allocate(left(nclass), right(nclass))
      left = 0.0_dp
      right = 0.0_dp
      do i = 1, size(indices)
         if ((ieee_is_nan(values(indices(i))) .and. .not. nan_right) .or. &
            (.not. ieee_is_nan(values(indices(i))) .and. values(indices(i)) <= cut)) then
            left(y(indices(i))) = left(y(indices(i))) + weights(i)
         else
            right(y(indices(i))) = right(y(indices(i))) + weights(i)
         end if
      end do
      ok = buckets_ok(left, right, min_bucket_class)
      if (.not. ok) then
         score = -huge(1.0_dp)
      else if (options%split_rule == RANGER_SPLIT_HELLINGER) then
         score = hellinger_score(left, right)
      else
         score = gini_score(left, class_weights) + gini_score(right, class_weights)
      end if
   end subroutine class_split_score

   subroutine best_categorical_class(values, y, indices, weights, ncat, nclass, class_weights, options, &
      min_bucket_class, rng, best_score, best_mask, best_nan_right, found)
      real(dp), intent(in) :: values(:), weights(:), class_weights(:)
      integer, intent(in) :: y(:), indices(:), ncat, nclass, min_bucket_class(:)
      type(ranger_options), intent(in) :: options
      type(ranger_rng_state), intent(inout) :: rng
      real(dp), intent(out) :: best_score
      logical, intent(out) :: best_mask(:), best_nan_right, found
      logical, allocatable :: mask(:)
      integer, allocatable :: present(:)
      integer :: r, ntry, npresent
      integer(i64) :: split_id, max_split_id
      real(dp) :: score_l, score_r
      logical :: ok_l, ok_r

      allocate(mask(ncat))
      best_score = -huge(1.0_dp)
      best_mask = .false.
      best_nan_right = .false.
      found = .false.

      if (options%split_rule == RANGER_SPLIT_EXTRATREES) then
         ntry = max(1, options%num_random_splits)
         do r = 1, ntry
            call random_extratrees_factor_mask(rng, values, indices, ncat, mask)
            call categorical_class_score(values, y, indices, weights, mask, .false., nclass, class_weights, &
               options, min_bucket_class, score_l, ok_l)
            call categorical_class_score(values, y, indices, weights, mask, .true., nclass, class_weights, &
               options, min_bucket_class, score_r, ok_r)
            call choose_mask_direction(score_l, ok_l, score_r, ok_r, mask, best_score, best_mask, &
               best_nan_right, found)
         end do
      else
         call collect_present_categories(values, indices, ncat, present, npresent)
         if (npresent < 2) return
         max_split_id = shiftl(1_i64, npresent - 1) - 1_i64
         do split_id = 1_i64, max_split_id
            call partition_mask_from_id(present, npresent, split_id, mask)
            call categorical_class_score(values, y, indices, weights, mask, .false., nclass, class_weights, &
               options, min_bucket_class, score_l, ok_l)
            call categorical_class_score(values, y, indices, weights, mask, .true., nclass, class_weights, &
               options, min_bucket_class, score_r, ok_r)
            call choose_mask_direction(score_l, ok_l, score_r, ok_r, mask, best_score, best_mask, &
               best_nan_right, found)
         end do
      end if
   end subroutine best_categorical_class

   subroutine categorical_class_score(values, y, indices, weights, mask, nan_right, nclass, class_weights, options, &
      min_bucket_class, score, ok)
      real(dp), intent(in) :: values(:), weights(:), class_weights(:)
      integer, intent(in) :: y(:), indices(:), nclass, min_bucket_class(:)
      logical, intent(in) :: mask(:), nan_right
      type(ranger_options), intent(in) :: options
      real(dp), intent(out) :: score
      logical, intent(out) :: ok
      real(dp), allocatable :: left(:), right(:)
      integer :: i, cat

      allocate(left(nclass), right(nclass))
      left = 0.0_dp
      right = 0.0_dp
      do i = 1, size(indices)
         if (ieee_is_nan(values(indices(i)))) then
            if (nan_right) then
               right(y(indices(i))) = right(y(indices(i))) + weights(i)
            else
               left(y(indices(i))) = left(y(indices(i))) + weights(i)
            end if
         else
            cat = nint(values(indices(i)))
            if (cat >= 1 .and. cat <= size(mask)) then
               if (mask(cat)) then
                  left(y(indices(i))) = left(y(indices(i))) + weights(i)
               else
                  right(y(indices(i))) = right(y(indices(i))) + weights(i)
               end if
            else
               right(y(indices(i))) = right(y(indices(i))) + weights(i)
            end if
         end if
      end do
      ok = buckets_ok(left, right, min_bucket_class)
      if (.not. ok) then
         score = -huge(1.0_dp)
      else if (options%split_rule == RANGER_SPLIT_HELLINGER) then
         score = hellinger_score(left, right)
      else
         score = gini_score(left, class_weights) + gini_score(right, class_weights)
      end if
   end subroutine categorical_class_score


   pure real(dp) function gini_score(counts, class_weights) result(score)
      real(dp), intent(in) :: counts(:), class_weights(:)
      real(dp) :: n
      n = sum(counts)
      if (n <= 0.0_dp) then
         score = 0.0_dp
      else
         score = sum(class_weights * counts * counts) / n
      end if
   end function gini_score

   pure real(dp) function hellinger_score(left, right) result(score)
      real(dp), intent(in) :: left(:), right(:)
      real(dp) :: total0, total1, tpr, fpr, a1, a2

      if (size(left) /= 2) then
         score = -huge(1.0_dp)
         return
      end if
      total0 = left(1) + right(1)
      total1 = left(2) + right(2)
      if (total0 <= 0.0_dp .or. total1 <= 0.0_dp) then
         score = -huge(1.0_dp)
         return
      end if
      tpr = right(2) / total1
      fpr = right(1) / total0
      a1 = sqrt(max(tpr, 0.0_dp)) - sqrt(max(fpr, 0.0_dp))
      a2 = sqrt(max(1.0_dp - tpr, 0.0_dp)) - sqrt(max(1.0_dp - fpr, 0.0_dp))
      score = sqrt(a1 * a1 + a2 * a2)
   end function hellinger_score

   logical function buckets_ok(left, right, min_bucket_class) result(ok)
      real(dp), intent(in) :: left(:), right(:)
      integer, intent(in) :: min_bucket_class(:)
      integer :: k

      ok = .true.
      if (sum(left) <= 0.0_dp .or. sum(right) <= 0.0_dp) then
         ok = .false.
         return
      end if
      if (size(min_bucket_class) == 1) then
         ok = sum(left) >= real(min_bucket_class(1), dp) .and. sum(right) >= real(min_bucket_class(1), dp)
      else
         do k = 1, min(size(left), size(min_bucket_class))
            if (left(k) < real(min_bucket_class(k), dp) .or. right(k) < real(min_bucket_class(k), dp)) then
               ok = .false.
               return
            end if
         end do
      end if
   end function buckets_ok

   subroutine choose_nan_direction(score_l, ok_l, score_r, ok_r, cut, best_score, best_cut, best_nan_right, found)
      real(dp), intent(in) :: score_l, score_r, cut
      logical, intent(in) :: ok_l, ok_r
      real(dp), intent(inout) :: best_score, best_cut
      logical, intent(inout) :: best_nan_right, found

      if (ok_l .and. (.not. found .or. score_l > best_score)) then
         found = .true.
         best_score = score_l
         best_cut = cut
         best_nan_right = .false.
      end if
      if (ok_r .and. (.not. found .or. score_r > best_score)) then
         found = .true.
         best_score = score_r
         best_cut = cut
         best_nan_right = .true.
      end if
   end subroutine choose_nan_direction

   subroutine choose_mask_direction(score_l, ok_l, score_r, ok_r, mask, best_score, best_mask, &
      best_nan_right, found)
      real(dp), intent(in) :: score_l, score_r
      logical, intent(in) :: ok_l, ok_r, mask(:)
      real(dp), intent(inout) :: best_score
      logical, intent(inout) :: best_mask(:), best_nan_right, found

      if (ok_l .and. (.not. found .or. score_l > best_score)) then
         found = .true.
         best_score = score_l
         best_mask = .false.
         best_mask(1:size(mask)) = mask
         best_nan_right = .false.
      end if
      if (ok_r .and. (.not. found .or. score_r > best_score)) then
         found = .true.
         best_score = score_r
         best_mask = .false.
         best_mask(1:size(mask)) = mask
         best_nan_right = .true.
      end if
   end subroutine choose_mask_direction

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

   subroutine predict_classification_tree(tree, x, ncat, prediction, probability, terminal_node)
      type(ranger_tree), intent(in) :: tree
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: ncat(:)
      integer, intent(out), optional :: prediction(:), terminal_node(:)
      real(dp), intent(out), optional :: probability(:,:)
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
         if (present(prediction)) prediction(i) = tree%node_class(node)
         if (present(probability)) probability(i, :) = tree%class_prob(:, node)
         if (present(terminal_node)) terminal_node(i) = node
      end do
   end subroutine predict_classification_tree

end module ranger_tree_classification
