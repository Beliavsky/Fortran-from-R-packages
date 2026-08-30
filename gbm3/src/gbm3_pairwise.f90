! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of gbm3 3.0.3; see NOTICE.md and upstream/.
module gbm3_pairwise
   use gbm3_kinds, only : dp
   use gbm3_constants
   use gbm3_math, only : argsort_real, argsort_desc_real, sigmoid, quiet_nan
   use gbm3_types, only : gbm_options, gbm_tree
   implicit none
   private
   public :: pairwise_working_response, pairwise_deviance, pairwise_fit_best_constant
   public :: pairwise_bag_improvement, pairwise_make_bag

contains

   subroutine pairwise_make_bag(group, bag_fraction, in_bag)
      integer, intent(in) :: group(:)
      real(dp), intent(in) :: bag_fraction
      logical, intent(out) :: in_bag(size(group))
      integer :: n_groups, target, seen, chosen, i, j
      real(dp) :: u

      in_bag = .false.
      if (size(group) == 0) return
      n_groups = count_groups(group)
      target = int(bag_fraction * real(n_groups, dp))
      if (target <= 0) target = 1
      seen = 0
      chosen = 0
      i = 1
      do while (i <= size(group))
         j = i + 1
         do while (j <= size(group))
            if (group(j) /= group(i)) exit
            j = j + 1
         end do
         if (chosen >= target) exit
         call random_number(u)
         if (u * real(n_groups - seen, dp) < real(target - chosen, dp)) then
            in_bag(i:j - 1) = .true.
            chosen = chosen + 1
         end if
         seen = seen + 1
         i = j
      end do
   end subroutine pairwise_make_bag

   integer function count_groups(group) result(n)
      integer, intent(in) :: group(:)
      integer :: i
      if (size(group) == 0) then
         n = 0
         return
      end if
      n = 1
      do i = 2, size(group)
         if (group(i) /= group(i - 1)) n = n + 1
      end do
   end function count_groups

   subroutine pairwise_working_response(y, group, offset, f, weight, in_bag, options, residual, hessian)
      real(dp), intent(in) :: y(:), offset(:), f(:), weight(:)
      integer, intent(in) :: group(:)
      logical, intent(in) :: in_bag(:)
      type(gbm_options), intent(in) :: options
      real(dp), intent(out) :: residual(size(y)), hessian(size(y))
      integer :: i, j
      real(dp), allocatable :: scores(:)

      residual = 0.0_dp
      hessian = 0.0_dp
      i = 1
      do while (i <= size(y))
         j = i + 1
         do while (j <= size(y))
            if (group(j) /= group(i)) exit
            j = j + 1
         end do
         if (in_bag(i)) then
            allocate(scores(j - i))
            scores = f(i:j - 1) + offset(i:j - 1)
            call compute_lambdas(group(i), y(i:j - 1), scores, weight(i:j - 1), options, &
                                 residual(i:j - 1), hessian(i:j - 1))
            deallocate(scores)
         end if
         i = j
      end do
   end subroutine pairwise_working_response

   subroutine compute_lambdas(group_id, y, scores, weight, options, residual, deriv)
      integer, intent(in) :: group_id
      real(dp), intent(in) :: y(:), scores(:), weight(:)
      type(gbm_options), intent(in) :: options
      real(dp), intent(inout) :: residual(:), deriv(:)
      integer, allocatable :: rank(:), item_at_rank(:)
      integer :: i, j, label_start, pairs
      real(dp) :: max_score, label_current, swap_cost, rho, lambda, d2, qnorm

      if (size(y) < 2 .or. weight(1) <= 0.0_dp) return
      max_score = metric_max(group_id, y, options)
      if (max_score <= 0.0_dp) return
      allocate(rank(size(y)), item_at_rank(size(y)))
      call rank_scores(scores, rank, item_at_rank)

      label_current = y(1)
      label_start = 0
      pairs = 0
      do j = 2, size(y)
         if (y(j) < label_current .or. y(j) > label_current) then
            label_start = j - 1
            label_current = y(j)
         end if
         do i = 1, label_start
            swap_cost = abs(metric_swap_cost(i, j, y, rank, item_at_rank, options))
            if (swap_cost > 0.0_dp) then
               pairs = pairs + 1
               rho = sigmoid(scores(j) - scores(i))
               lambda = swap_cost * rho
               residual(i) = residual(i) + lambda
               residual(j) = residual(j) - lambda
               d2 = lambda * (1.0_dp - rho)
               deriv(i) = deriv(i) + d2
               deriv(j) = deriv(j) + d2
            end if
         end do
      end do
      if (pairs > 0) then
         qnorm = 1.0_dp / (max_score * real(pairs, dp))
         residual = residual * qnorm
         deriv = deriv * qnorm
      end if
   end subroutine compute_lambdas

   subroutine rank_scores(scores, rank, item_at_rank)
      real(dp), intent(in) :: scores(:)
      integer, intent(out) :: rank(size(scores)), item_at_rank(size(scores))
      real(dp), allocatable :: noisy(:)
      integer, allocatable :: ord(:)
      real(dp) :: u
      integer :: i
      allocate(noisy(size(scores)), ord(size(scores)))
      do i = 1, size(scores)
         call random_number(u)
         noisy(i) = scores(i) + 1.0e-10_dp * (u - 0.5_dp)
      end do
      call argsort_desc_real(noisy, ord)
      do i = 1, size(scores)
         item_at_rank(i) = ord(i)
         rank(ord(i)) = i
      end do
   end subroutine rank_scores

   pure logical function any_pairs(y) result(ok)
      real(dp), intent(in) :: y(:)
      ok = size(y) >= 2
      if (.not. ok) return
      ok = y(1) > 0.0_dp .and. (y(size(y)) < y(1) .or. y(size(y)) > y(1))
   end function any_pairs

   real(dp) function metric_max(group_id, y, options) result(v)
      integer, intent(in) :: group_id
      real(dp), intent(in) :: y(:)
      type(gbm_options), intent(in) :: options
      integer :: i, cutoff
      real(dp) :: rw
      associate(dummy => group_id)
      end associate
      if (.not. any_pairs(y)) then
         v = 0.0_dp
         return
      end if
      select case (options%pairwise_metric)
      case (GBM_METRIC_CONC)
         v = real(pair_count(y), dp)
      case (GBM_METRIC_NDCG)
         cutoff = options%pairwise_max_rank
         if (cutoff <= 0) cutoff = size(y)
         v = 0.0_dp
         do i = 1, min(size(y), cutoff)
            if (y(i) <= 0.0_dp) exit
            rw = log(2.0_dp) / log(real(i + 1, dp))
            v = v + y(i) * rw
         end do
      case (GBM_METRIC_MAP, GBM_METRIC_MRR)
         v = 1.0_dp
      case default
         v = 0.0_dp
      end select
   end function metric_max

   integer function pair_count(y) result(pairs)
      real(dp), intent(in) :: y(:)
      real(dp) :: label_current
      integer :: j, label_end
      if (.not. any_pairs(y)) then
         pairs = 0
         return
      end if
      label_current = y(1)
      label_end = 0
      pairs = 0
      do j = 2, size(y)
         if (y(j) < label_current .or. y(j) > label_current) then
            label_end = j - 1
            label_current = y(j)
         end if
         pairs = pairs + label_end
      end do
   end function pair_count

   real(dp) function metric_measure(y, rank, item_at_rank, options) result(v)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: rank(:), item_at_rank(:)
      type(gbm_options), intent(in) :: options
      integer :: i, j, label_end, good, cutoff, pos, top
      integer, allocatable :: posranks(:)
      real(dp) :: label_current
      associate(dummy => item_at_rank)
      end associate
      select case (options%pairwise_metric)
      case (GBM_METRIC_CONC)
         label_current = y(1)
         label_end = 0
         good = 0
         do j = 2, size(y)
            if (y(j) < label_current .or. y(j) > label_current) then
               label_end = j - 1
               label_current = y(j)
            end if
            do i = 1, label_end
               if (rank(i) < rank(j)) good = good + 1
            end do
         end do
         v = real(good, dp)
      case (GBM_METRIC_NDCG)
         cutoff = options%pairwise_max_rank
         if (cutoff <= 0) cutoff = size(y)
         v = 0.0_dp
         do i = 1, size(y)
            if (rank(i) <= cutoff) v = v + y(i) * log(2.0_dp) / log(real(rank(i) + 1, dp))
         end do
      case (GBM_METRIC_MRR)
         cutoff = options%pairwise_max_rank
         if (cutoff <= 0) cutoff = size(y)
         top = size(y) + 1
         do i = 1, size(y)
            if (y(i) <= 0.0_dp) exit
            top = min(top, rank(i))
         end do
         if (top > cutoff) then
            v = 0.0_dp
         else
            v = 1.0_dp / real(top, dp)
         end if
      case (GBM_METRIC_MAP)
         pos = 0
         do i = 1, size(y)
            if (y(i) <= 0.0_dp) exit
            pos = pos + 1
         end do
         if (pos == 0) then
            v = 0.0_dp
         else
            allocate(posranks(pos))
            posranks = rank(1:pos)
            call sort_int(posranks)
            v = 0.0_dp
            do i = 1, pos
               v = v + real(i, dp) / real(posranks(i), dp)
            end do
            v = v / real(pos, dp)
         end if
      case default
         v = 0.0_dp
      end select
   end function metric_measure

   real(dp) function metric_swap_cost(item_better, item_worse, y, rank, item_at_rank, options) result(v)
      integer, intent(in) :: item_better, item_worse
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: rank(:), item_at_rank(:)
      type(gbm_options), intent(in) :: options
      integer :: rb, rw, rank_upper, rank_lower, r, diff
      integer :: pos, top, cutoff, npos_pos, npos_neg, lo, hi, j
      integer, allocatable :: posranks(:)
      real(dp) :: resp_upper, resp_lower, score_diff, cur, negm, before, after, sg

      select case (options%pairwise_metric)
      case (GBM_METRIC_NDCG)
         v = (rank_weight(rank(item_better), options) - rank_weight(rank(item_worse), options)) * &
             (y(item_better) - y(item_worse))
      case (GBM_METRIC_CONC)
         rb = rank(item_better)
         rw = rank(item_worse)
         if (rb > rw) then
            rank_upper = rw
            rank_lower = rb
            resp_upper = y(item_worse)
            resp_lower = y(item_better)
            diff = 1
         else
            rank_upper = rb
            rank_lower = rw
            resp_upper = y(item_better)
            resp_lower = y(item_worse)
            diff = -1
         end if
         do r = rank_upper + 1, rank_lower - 1
            score_diff = y(item_at_rank(r)) - resp_lower
            if (score_diff < 0.0_dp) then
               diff = diff + 1
            else if (score_diff > 0.0_dp) then
               diff = diff - 1
            end if
            score_diff = y(item_at_rank(r)) - resp_upper
            if (score_diff < 0.0_dp) then
               diff = diff - 1
            else if (score_diff > 0.0_dp) then
               diff = diff + 1
            end if
         end do
         v = real(diff, dp)
      case (GBM_METRIC_MRR)
         pos = positive_count(y)
         if (pos == 0 .or. pos >= size(y)) then
            v = 0.0_dp
            return
         end if
         top = minval(rank(1:pos))
         rb = rank(item_better)
         rw = rank(item_worse)
         cutoff = options%pairwise_max_rank
         if (cutoff <= 0) cutoff = size(y)
         if (top > cutoff) then
            cur = 0.0_dp
         else
            cur = 1.0_dp / real(top, dp)
         end if
         if (rw > cutoff) then
            negm = 0.0_dp
         else
            negm = 1.0_dp / real(rw, dp)
         end if
         if (rw < top .or. rb == top) then
            v = negm - cur
         else
            v = 0.0_dp
         end if
      case (GBM_METRIC_MAP)
         pos = positive_count(y)
         if (pos == 0) then
            v = 0.0_dp
            return
         end if
         allocate(posranks(pos))
         posranks = rank(1:pos)
         call sort_int(posranks)
         rb = rank(item_better)
         rw = rank(item_worse)
         npos_pos = count(posranks <= rb)
         npos_neg = count(posranks <= rw)
         before = real(npos_pos, dp) / real(rb, dp)
         if (rw > rb) then
            sg = -1.0_dp
            lo = npos_pos + 1
            hi = npos_neg - 1
            after = real(npos_neg, dp) / real(rw, dp)
         else
            sg = 1.0_dp
            lo = npos_neg + 1
            hi = npos_pos - 1
            after = real(npos_neg + 1, dp) / real(rw, dp)
         end if
         v = after - before
         do j = lo, hi
            if (j >= 1 .and. j <= pos) v = v + sg / real(posranks(j), dp)
         end do
         v = v / real(pos, dp)
      case default
         v = 0.0_dp
      end select
   end function metric_swap_cost

   pure real(dp) function rank_weight(rank, options) result(w)
      integer, intent(in) :: rank
      type(gbm_options), intent(in) :: options
      integer :: cutoff
      cutoff = options%pairwise_max_rank
      if (cutoff <= 0) cutoff = huge(0)
      if (rank > cutoff) then
         w = 0.0_dp
      else
         w = log(2.0_dp) / log(real(rank + 1, dp))
      end if
   end function rank_weight

   integer function positive_count(y) result(n)
      real(dp), intent(in) :: y(:)
      integer :: i
      n = 0
      do i = 1, size(y)
         if (y(i) <= 0.0_dp) exit
         n = n + 1
      end do
   end function positive_count

   subroutine sort_int(x)
      integer, intent(inout) :: x(:)
      integer :: i, j, key
      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine sort_int

   real(dp) function pairwise_deviance(y, group, offset, f, weight, options) result(dev)
      real(dp), intent(in) :: y(:), offset(:), f(:), weight(:)
      integer, intent(in) :: group(:)
      type(gbm_options), intent(in) :: options
      integer :: i, j
      integer, allocatable :: rank(:), item_at_rank(:)
      real(dp), allocatable :: scores(:)
      real(dp) :: loss, sw, mx

      if (size(y) == 0) then
         dev = 0.0_dp
         return
      end if
      loss = 0.0_dp
      sw = 0.0_dp
      i = 1
      do while (i <= size(y))
         j = i + 1
         do while (j <= size(y))
            if (group(j) /= group(i)) exit
            j = j + 1
         end do
         mx = metric_max(group(i), y(i:j - 1), options)
         if (mx > 0.0_dp) then
            allocate(scores(j - i), rank(j - i), item_at_rank(j - i))
            scores = f(i:j - 1) + offset(i:j - 1)
            call rank_scores(scores, rank, item_at_rank)
            loss = loss + weight(i) * metric_measure(y(i:j - 1), rank, item_at_rank, options) / mx
            sw = sw + weight(i)
            deallocate(scores, rank, item_at_rank)
         end if
         i = j
      end do
      if (sw <= 0.0_dp) then
         dev = quiet_nan()
      else
         dev = 1.0_dp - loss / sw
      end if
   end function pairwise_deviance

   subroutine pairwise_fit_best_constant(weight, residual, hessian, in_bag, assignment, tree)
      real(dp), intent(in) :: weight(:), residual(:), hessian(:)
      logical, intent(in) :: in_bag(:)
      integer, intent(in) :: assignment(:)
      type(gbm_tree), intent(inout) :: tree
      real(dp), allocatable :: num(:), den(:)
      integer :: i, node
      allocate(num(tree%n_nodes), den(tree%n_nodes))
      num = 0.0_dp
      den = 0.0_dp
      do i = 1, size(weight)
         if (.not. in_bag(i)) cycle
         node = assignment(i)
         num(node) = num(node) + weight(i) * residual(i)
         den(node) = den(node) + weight(i) * hessian(i)
      end do
      do node = 1, tree%n_nodes
         if (.not. tree%nodes(node)%is_terminal) cycle
         if (den(node) <= 0.0_dp) then
            tree%nodes(node)%prediction = 0.0_dp
         else
            tree%nodes(node)%prediction = num(node) / den(node)
         end if
      end do
   end subroutine pairwise_fit_best_constant

   real(dp) function pairwise_bag_improvement(y, group, offset, f, weight, delta, in_bag, options) result(improvement)
      real(dp), intent(in) :: y(:), offset(:), f(:), weight(:), delta(:)
      integer, intent(in) :: group(:)
      logical, intent(in) :: in_bag(:)
      type(gbm_options), intent(in) :: options
      integer :: i, j
      integer, allocatable :: rank(:), item_at_rank(:)
      real(dp), allocatable :: oldscores(:), newscores(:)
      real(dp) :: mx, oldm, newm, val, sw

      val = 0.0_dp
      sw = 0.0_dp
      i = 1
      do while (i <= size(y))
         j = i + 1
         do while (j <= size(y))
            if (group(j) /= group(i)) exit
            j = j + 1
         end do
         if (.not. in_bag(i)) then
            mx = metric_max(group(i), y(i:j - 1), options)
            if (mx > 0.0_dp) then
               allocate(oldscores(j - i), newscores(j - i), rank(j - i), item_at_rank(j - i))
               oldscores = f(i:j - 1) + offset(i:j - 1)
               newscores = oldscores + options%shrinkage * delta(i:j - 1)
               call rank_scores(oldscores, rank, item_at_rank)
               oldm = metric_measure(y(i:j - 1), rank, item_at_rank, options)
               call rank_scores(newscores, rank, item_at_rank)
               newm = metric_measure(y(i:j - 1), rank, item_at_rank, options)
               val = val + weight(i) * (newm - oldm) / mx
               sw = sw + weight(i)
               deallocate(oldscores, newscores, rank, item_at_rank)
            end if
         end if
         i = j
      end do
      if (sw <= 0.0_dp) then
         improvement = quiet_nan()
      else
         improvement = val / sw
      end if
   end function pairwise_bag_improvement

end module gbm3_pairwise
