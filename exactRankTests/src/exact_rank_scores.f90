! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation of exactRankTests/R/cscores.R and src/irank.c.
module exact_rank_scores
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use r_kinds, only : dp
   use r_sorting, only : r_average_ranks
   use r_ordering, only : r_order
   use r_distributions, only : r_qnorm
   use exact_rank_distribution, only : rank_ok, rank_invalid, rank_numeric
   implicit none
   private
   public :: cscores, cscores_factor, cscores_surv, irank
contains
   pure subroutine irank(x, ranks, status, order)
      real(dp), intent(in) :: x(:) !! Finite observations, with exact equality defining ties.
      integer, allocatable, intent(out) :: ranks(:) !! Maximum ranks: count of observations <= each observation.
      integer, intent(out) :: status !! rank_ok on success; rank_invalid for nonfinite data or bad order.
      integer, intent(in), optional :: order(:) !! Optional one-based permutation sorting x in nondecreasing order.
      integer, allocatable :: indices(:), seen(:)
      integer :: i, j, n
      status = rank_invalid
      if (any(.not. ieee_is_finite(x))) return
      n = size(x)
      if (present(order)) then
         if (size(order) /= n) return
         if (any(order < 1) .or. any(order > n)) return
         allocate(seen(n))
         seen = 0
         do i = 1, n
            seen(order(i)) = seen(order(i)) + 1
         end do
         if (any(seen /= 1)) return
         indices = order
         do i = 2, n
            if (x(indices(i)) < x(indices(i - 1))) return
         end do
      else
         call r_order(x, indices)
      end if
      allocate(ranks(n))
      i = 1
      do while (i <= n)
         j = i
         do while (j < n)
            if (x(indices(j + 1)) /= x(indices(i))) exit
            j = j + 1
         end do
         ranks(indices(i:j)) = j
         i = j + 1
      end do
      status = rank_ok
   end subroutine irank

   pure subroutine cscores(y, scores, status, score_type, integer_scores, max_score)
      real(dp), intent(in) :: y(:) !! Finite numerical observations; an empty vector returns an empty score vector.
      real(dp), allocatable, intent(out) :: scores(:) !! Scores in observation order, length size(y).
      integer, intent(out) :: status !! rank_ok on success; invalid type or nonfinite data yields rank_invalid.
      character(*), intent(in), optional :: score_type !! Data (default), Wilcoxon, NormalQuantile, AnsariBradley, Median, Savage, ConSal.
      logical, intent(in), optional :: integer_scores !! Apply upstream integer rescaling; default false.
      integer, intent(in), optional :: max_score !! Positive rescaling bound, default size(y), used by non-midrank scores.
      real(dp), allocatable :: ranks(:)
      integer, allocatable :: event(:)
      character(16) :: kind
      integer :: n
      status = rank_invalid
      n = size(y)
      if (any(.not. ieee_is_finite(y))) return
      kind = 'Data'
      if (present(score_type)) kind = score_type
      if (kind == 'Data') then
         scores = y
      else
         call r_average_ranks(y, ranks)
         select case (kind)
         case ('Wilcoxon')
            scores = ranks
         case ('NormalQuantile')
            scores = r_qnorm(ranks / real(n + 1, dp))
         case ('AnsariBradley', 'Ansari')
            scores = min(ranks, real(n + 1, dp) - ranks)
         case ('Median')
            scores = merge(1.0_dp, 0.0_dp, ranks > real(n + 1, dp) / 2.0_dp)
         case ('Savage')
            allocate(event(n))
            event = 1
            call cscores_surv(y, event, scores, status)
            if (status /= rank_ok) return
         case ('ConSal')
            scores = (ranks / real(n + 1, dp)) ** 4
         case default
            return
         end select
      end if
      call integer_rescale(scores, status, integer_scores, max_score, .true.)
   end subroutine cscores

   pure subroutine cscores_factor(levels, scores, status)
      integer, intent(in) :: levels(:) !! Factor codes 1 or 2; R labels and level construction are caller responsibilities.
      real(dp), allocatable, intent(out) :: scores(:) !! Numerical copies of the one-based factor codes.
      integer, intent(out) :: status !! rank_ok if every level code is 1 or 2, rank_invalid otherwise.
      status = rank_invalid
      if (any(levels < 1) .or. any(levels > 2)) return
      scores = real(levels, dp)
      status = rank_ok
   end subroutine cscores_factor

   pure subroutine cscores_surv(time, event, scores, status, integer_scores, max_score)
      real(dp), intent(in) :: time(:) !! Finite observed event/censoring times.
      integer, intent(in) :: event(:) !! Event indicators, same length as time, each 0 (censored) or 1 (event).
      real(dp), allocatable, intent(out) :: scores(:) !! Log-rank scores in the original observation order.
      integer, intent(out) :: status !! rank_ok on success; rank_invalid for malformed input.
      logical, intent(in), optional :: integer_scores !! Round scaled scores to integers, default false.
      integer, intent(in), optional :: max_score !! Positive rescaling bound, default size(time).
      integer, allocatable :: ranks(:), order(:)
      real(dp), allocatable :: cumulative(:)
      integer :: i, n
      status = rank_invalid
      n = size(time)
      if (size(event) /= n .or. any(event < 0) .or. any(event > 1)) return
      call irank(time, ranks, status)
      if (status /= rank_ok) return
      call r_order(time, order)
      allocate(cumulative(n), scores(n))
      do i = 1, n
         cumulative(i) = real(event(order(i)), dp) / real(n - ranks(order(i)) + 1, dp)
         if (i > 1) cumulative(i) = cumulative(i) + cumulative(i - 1)
      end do
      scores = real(event, dp) - cumulative(ranks)
      call integer_rescale(scores, status, integer_scores, max_score, .false.)
   end subroutine cscores_surv

   pure subroutine integer_rescale(scores, status, integer_scores, max_score, midrank_shortcut)
      real(dp), intent(inout) :: scores(:) !! Finite scores to rescale in place when requested.
      integer, intent(out) :: status !! rank_ok on success; rank_numeric if rescaling has zero denominator.
      logical, intent(in), optional :: integer_scores !! Whether rescaling is requested, default false.
      integer, intent(in), optional :: max_score !! Positive requested maximum, default number of scores.
      logical, intent(in) :: midrank_shortcut !! Preserve cscores.default's doubling of integer and half-integer scores.
      integer :: bound, i
      real(dp) :: value, a
      logical :: transform
      status = rank_invalid
      bound = max(1, size(scores))
      if (present(max_score)) bound = max_score
      if (bound < 1) return
      transform = .false.
      if (present(integer_scores)) transform = integer_scores
      status = rank_ok
      if (.not. transform .or. size(scores) == 0) return
      if (maxval(abs(scores)) > 2.0_dp ** 50) then
         status = rank_numeric
         return
      end if
      if (midrank_shortcut) then
         if (all(2.0_dp * scores == anint(2.0_dp * scores))) then
            scores = 2.0_dp * scores
            return
         end if
      end if
      if (maxval(scores) == 0.0_dp) then
         status = rank_numeric
         return
      end if
      scores = scores * (real(bound, dp) / maxval(scores))
      if (any(.not. ieee_is_finite(scores)) .or. maxval(abs(scores)) > 2.0_dp ** 50) then
         status = rank_numeric
         return
      end if
      do i = 1, size(scores)
         value = abs(scores(i))
         a = aint(value)
         if (value - a > 0.5_dp .or. (value - a == 0.5_dp .and. modulo(a, 2.0_dp) /= 0.0_dp)) a = a + 1.0_dp
         scores(i) = sign(a, scores(i))
      end do
   end subroutine integer_rescale
end module exact_rank_scores
