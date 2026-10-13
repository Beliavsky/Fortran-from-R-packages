! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from exactRankTests by Torsten Hothorn and Kurt Hornik.
! Copyright (C) 2000-2003 Torsten Hothorn, upstream algorithms.
module exact_rank_distribution
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use r_ordering, only : r_order, r_sort_values_in_place
   implicit none
   private
   integer, parameter, public :: rank_ok = 0, rank_invalid = 1, rank_limit = 2, rank_allocation = 3
   integer, parameter, public :: rank_numeric = 4
   integer(int64), parameter :: default_cells = 5000000_int64
   type, public :: rank_distribution
      real(dp), allocatable :: support(:)
      real(dp), allocatable :: probability(:)
      real(dp) :: expectation = 0.0_dp
      real(dp) :: factor = 1.0_dp
      real(dp) :: rounding_error = 0.0_dp
      logical :: paired = .false.
      integer :: status = rank_invalid
   end type rank_distribution
   public :: dp, permutation_distribution, simulate_distribution, dperm, pperm, qperm, rperm
contains
   pure elemental real(dp) function round_even(x) result(value)
      real(dp), intent(in) :: x !! Finite real score within the exactly representable integer range.
      real(dp) :: a, f
      a = floor(abs(x), kind = int64)
      f = abs(x) - a
      if (f > 0.5_dp .or. (f == 0.5_dp .and. modulo(a, 2.0_dp) /= 0.0_dp)) a = a + 1.0_dp
      value = sign(a, x)
   end function round_even

   pure subroutine permutation_distribution(scores, m, distribution, status, paired, factor, tolerance, max_cells)
      real(dp), intent(in) :: scores(:) !! Finite scores; signed scores are allowed for independent samples.
      integer, intent(in) :: m !! Number sampled without replacement; 1 <= m <= size(scores).
      type(rank_distribution), intent(out) :: distribution !! Sorted attainable statistics and normalized probabilities.
      integer, intent(out) :: status !! rank_ok on success; rank_limit prevents an oversized dynamic-programming table.
      logical, intent(in), optional :: paired !! Sign-flip distribution only when m equals size(scores); defaults to that case.
      real(dp), intent(in), optional :: factor !! Positive lattice multiplier; scores are rounded with R's ties-to-even rule.
      real(dp), intent(in), optional :: tolerance !! Automatic lattice approximation error target, default 0.01.
      integer(int64), intent(in), optional :: max_cells !! Maximum table cells, default 5,000,000 (approximately 40 MB).
      integer(int64), allocatable :: integer_scores(:)
      real(dp), allocatable :: work(:, :), reduced(:), sorted(:), approximated(:)
      integer, allocatable :: ordering(:)
      integer(int64) :: base, divisor, span64, limit, cap, t, a, b, total64
      integer :: n, i, j, s, v, span, rows, top, count_support, k, astat, grid, trial
      real(dp) :: multiplier, target, best, error, criterion, lo, hi, total, shift, weight
      logical :: sign_flip
      status = rank_invalid
      n = size(scores)
      if (n < 1 .or. m < 1 .or. m > n) return
      if (any(.not. ieee_is_finite(scores))) return
      if (maxval(abs(scores)) > huge(1.0_dp) / real(n, dp)) then
         status = rank_numeric
         return
      end if
      sign_flip = n == m
      if (present(paired)) sign_flip = sign_flip .and. paired
      if (sign_flip .and. any(scores < 0.0_dp)) return
      target = 0.01_dp
      if (present(tolerance)) target = tolerance
      if (.not. ieee_is_finite(target) .or. target < 0.0_dp) return
      limit = default_cells
      if (present(max_cells)) limit = max_cells
      if (limit < 1_int64) return
      multiplier = 1.0_dp
      if (present(factor)) then
         multiplier = factor
         if (.not. ieee_is_finite(multiplier) .or. multiplier <= 0.0_dp) return
      else
         if (maxval(abs(scores)) > 2.0_dp ** 50) then
            status = rank_limit
            return
         end if
         if (any(scores /= round_even(scores))) then
            if (all(2.0_dp * scores == round_even(2.0_dp * scores))) then
               multiplier = 2.0_dp
            else
               ! The upstream search minimizes |rounding error - target|.
               reduced = scores - minval(scores) + 1.0_dp
               call r_order(reduced, ordering)
               sorted = reduced(ordering(n:1:-1))
               total = sum(sorted) - minval(sorted)
               if (.not. ieee_is_finite(total)) return
               total = max(total, 1.0_dp)
               grid = int(min(10000.0_dp, max(1.0_dp, 10000.0_dp / total)))
               best = huge(1.0_dp)
               do trial = 1, grid
                  error = sum(abs(real(trial, dp) * sorted(:m) - &
                     round_even(real(trial, dp) * sorted(:m)))) / real(trial, dp)
                  criterion = abs(error - target)
                  if (criterion < best) then
                     best = criterion
                     multiplier = real(trial, dp)
                  end if
               end do
            end if
         end if
      end if
      if (multiplier >= 1.0_dp) then
         if (maxval(abs(scores)) > (2.0_dp ** 50) / multiplier) then
            status = rank_limit
            return
         end if
      else
         if (maxval(abs(scores)) * multiplier > 2.0_dp ** 50) then
            status = rank_limit
            return
         end if
      end if
      approximated = round_even(scores * multiplier)
      distribution%rounding_error = sum(abs(scores - approximated / multiplier))
      integer_scores = int(approximated, int64)
      base = 0_int64
      if (.not. sign_flip) base = minval(integer_scores)
      integer_scores = integer_scores - base
      divisor = 0_int64
      do i = 1, n
         a = divisor
         b = integer_scores(i)
         do while (b /= 0_int64)
            t = modulo(a, b)
            a = b
            b = t
         end do
         divisor = a
      end do
      divisor = max(1_int64, divisor)
      integer_scores = integer_scores / divisor
      rows = m
      if (sign_flip) rows = 0
      cap = limit / (int(rows, int64) + 1_int64) - 1_int64
      if (cap < 0_int64) then
         status = rank_limit
         return
      end if
      reduced = real(integer_scores, dp)
      call r_sort_values_in_place(reduced)
      span64 = 0_int64
      top = m
      if (sign_flip) top = n
      do i = n - top + 1, n
         t = int(reduced(i), int64)
         if (t > cap - span64) then
            status = rank_limit
            return
         end if
         span64 = span64 + t
      end do
      if (span64 > int(huge(span) - 1, int64)) then
         status = rank_limit
         return
      end if
      span = int(span64)
      allocate(work(0:span, 0:rows), stat = astat)
      if (astat /= 0) then
         status = rank_allocation
         return
      end if
      work = 0.0_dp
      work(0, 0) = 1.0_dp
      total64 = 0_int64
      do i = 1, n
         v = int(min(integer_scores(i), span64 + 1_int64))
         if (sign_flip) then
            do s = min(span, int(total64) + v), 0, -1
               weight = 0.5_dp * work(s, 0)
               if (s >= v) weight = weight + 0.5_dp * work(s - v, 0)
               work(s, 0) = weight
            end do
            total64 = min(span64, total64 + int(v, int64))
         else
            do j = min(i, m), 1, -1
               lo = real(i - j, dp) / real(i, dp)
               hi = real(j, dp) / real(i, dp)
               do s = span, 0, -1
                  weight = lo * work(s, j)
                  if (s >= v) weight = weight + hi * work(s - v, j - 1)
                  work(s, j) = weight
               end do
            end do
         end if
      end do
      count_support = count(work(:, rows) > 0.0_dp)
      total = sum(work(:, rows))
      if (count_support < 1 .or. .not. ieee_is_finite(total) .or. total <= 0.0_dp) then
         status = rank_numeric
         return
      end if
      allocate(distribution%support(count_support), distribution%probability(count_support), stat = astat)
      if (astat /= 0) then
         status = rank_allocation
         return
      end if
      shift = real(base, dp) * real(m, dp)
      k = 0
      do s = 0, span
         if (work(s, rows) <= 0.0_dp) cycle
         k = k + 1
         distribution%support(k) = (real(s, dp) * real(divisor, dp) + shift) / multiplier
         distribution%probability(k) = work(s, rows) / total
      end do
      if (sign_flip) then
         distribution%expectation = sum(scores * 0.5_dp)
      else
         distribution%expectation = (real(m, dp) / real(n, dp)) * sum(scores)
      end if
      distribution%factor = multiplier
      distribution%paired = sign_flip
      distribution%status = rank_ok
      status = rank_ok
   end subroutine permutation_distribution

   pure elemental real(dp) function dperm(x, distribution) result(probability)
      real(dp), intent(in) :: x !! Statistic at which to evaluate mass; a non-support value has mass zero.
      type(rank_distribution), intent(in) :: distribution !! Successfully constructed exact or simulated distribution.
      integer :: i
      probability = ieee_value(0.0_dp, ieee_quiet_nan)
      if (distribution%status /= rank_ok .or. ieee_is_nan(x)) return
      probability = 0.0_dp
      do i = 1, size(distribution%support)
         if (distribution%support(i) == x) then
            probability = distribution%probability(i)
            return
         end if
      end do
   end function dperm

   pure elemental real(dp) function pperm(q, distribution, alternative) result(probability)
      real(dp), intent(in) :: q !! Observed statistic; one-sided tails include equality.
      type(rank_distribution), intent(in) :: distribution !! Successfully constructed exact or simulated distribution.
      character(*), intent(in), optional :: alternative !! 'less' (default), 'greater', or 'two.sided'.
      character(9) :: alt
      real(dp) :: lower, upper, distance
      probability = ieee_value(0.0_dp, ieee_quiet_nan)
      if (distribution%status /= rank_ok .or. ieee_is_nan(q)) return
      alt = 'less'
      if (present(alternative)) then
         if (len_trim(alternative) > len(alt)) return
         alt = alternative
      end if
      lower = sum(distribution%probability, mask = distribution%support <= q)
      upper = sum(distribution%probability, mask = distribution%support >= q)
      select case (alt)
      case ('less')
         probability = lower
      case ('greater')
         probability = upper
      case ('two.sided')
         if (distribution%paired) then
            probability = 2.0_dp * min(lower, upper)
         else
            distance = abs(q - distribution%expectation)
            probability = sum(distribution%probability, &
               mask = abs(distribution%support - distribution%expectation) >= distance)
         end if
      case default
         return
      end select
      probability = min(1.0_dp, probability)
   end function pperm

   pure elemental real(dp) function qperm(p, distribution) result(quantile)
      real(dp), intent(in) :: p !! Cumulative probability in [0,1]; endpoints return extreme support values.
      type(rank_distribution), intent(in) :: distribution !! Successfully constructed exact or simulated distribution.
      integer :: i
      real(dp) :: cumulative
      quantile = ieee_value(0.0_dp, ieee_quiet_nan)
      if (distribution%status /= rank_ok .or. .not. ieee_is_finite(p)) return
      if (p < 0.0_dp .or. p > 1.0_dp) return
      cumulative = 0.0_dp
      do i = 1, size(distribution%support)
         cumulative = cumulative + distribution%probability(i)
         if (cumulative >= p) then
            quantile = distribution%support(i)
            return
         end if
      end do
      quantile = distribution%support(size(distribution%support))
   end function qperm

   pure subroutine rperm(scores, m, uniforms, draws, status)
      real(dp), intent(in) :: scores(:) !! Finite population scores to sample without replacement.
      integer, intent(in) :: m !! Sample size, between 1 and the population size.
      real(dp), intent(in) :: uniforms(:, :) !! Shape (m, number of draws), independent uniforms in [0,1).
      real(dp), allocatable, intent(out) :: draws(:) !! One sampled score sum per column of uniforms.
      integer, intent(out) :: status !! rank_ok on success; rank_invalid for bad dimensions or values.
      integer, allocatable :: indices(:)
      integer :: i, j, k, swap, n
      status = rank_invalid
      n = size(scores)
      if (m < 1 .or. m > n .or. size(uniforms, 1) /= m) return
      if (any(.not. ieee_is_finite(scores)) .or. any(.not. ieee_is_finite(uniforms))) return
      if (any(uniforms < 0.0_dp) .or. any(uniforms >= 1.0_dp)) return
      if (maxval(abs(scores)) > huge(1.0_dp) / real(m, dp)) then
         status = rank_numeric
         return
      end if
      allocate(draws(size(uniforms, 2)), indices(n))
      do j = 1, size(draws)
         indices = [(i, i = 1, n)]
         draws(j) = 0.0_dp
         do i = 1, m
            k = i + int(uniforms(i, j) * real(n - i + 1, dp))
            swap = indices(i)
            indices(i) = indices(k)
            indices(k) = swap
            draws(j) = draws(j) + scores(indices(i))
         end do
      end do
      if (any(.not. ieee_is_finite(draws))) then
         status = rank_numeric
         return
      end if
      status = rank_ok
   end subroutine rperm

   pure subroutine simulate_distribution(scores, m, uniforms, distribution, status, paired)
      real(dp), intent(in) :: scores(:) !! Finite scores; sign-flip scores must be nonnegative.
      integer, intent(in) :: m !! Sample size between 1 and size(scores).
      real(dp), intent(in) :: uniforms(:, :) !! Shape (m, B); independent [0,1) uniforms, B >= 1.
      type(rank_distribution), intent(out) :: distribution !! Empirical mass distribution of simulated statistics.
      integer, intent(out) :: status !! rank_ok on success; nonzero on invalid or overflowing input.
      logical, intent(in), optional :: paired !! Default size(scores)==m; false gives fixed-size sampling even for m==N.
      real(dp), allocatable :: draws(:)
      integer :: i, j, k, n, b
      logical :: sign_flip
      status = rank_invalid
      n = size(scores)
      b = size(uniforms, 2)
      if (m < 1 .or. m > n .or. size(uniforms, 1) /= m .or. b < 1) return
      if (any(.not. ieee_is_finite(scores)) .or. any(.not. ieee_is_finite(uniforms))) return
      if (any(uniforms < 0.0_dp) .or. any(uniforms >= 1.0_dp)) return
      if (maxval(abs(scores)) > huge(1.0_dp) / real(n, dp)) then
         status = rank_numeric
         return
      end if
      sign_flip = m == n
      if (present(paired)) sign_flip = sign_flip .and. paired
      if (sign_flip) then
         if (any(scores < 0.0_dp)) return
         allocate(draws(b))
         do i = 1, b
            draws(i) = sum(scores, mask = uniforms(:, i) < 0.5_dp)
         end do
      else
         call rperm(scores, m, uniforms, draws, status)
         if (status /= rank_ok) return
      end if
      if (any(.not. ieee_is_finite(draws))) then
         status = rank_numeric
         return
      end if
      call r_sort_values_in_place(draws)
      k = 1
      do i = 2, b
         if (draws(i) /= draws(i - 1)) k = k + 1
      end do
      allocate(distribution%support(k), distribution%probability(k))
      k = 0
      i = 1
      do while (i <= b)
         j = i
         do while (j < b)
            if (draws(j + 1) /= draws(i)) exit
            j = j + 1
         end do
         k = k + 1
         distribution%support(k) = draws(i)
         distribution%probability(k) = real(j - i + 1, dp) / real(b, dp)
         i = j + 1
      end do
      distribution%expectation = (real(m, dp) / real(n, dp)) * sum(scores)
      if (sign_flip) distribution%expectation = sum(scores * 0.5_dp)
      distribution%paired = sign_flip
      distribution%status = rank_ok
      status = rank_ok
   end subroutine simulate_distribution
end module exact_rank_distribution
