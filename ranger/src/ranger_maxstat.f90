! SPDX-License-Identifier: GPL-3.0-only
! Modern Fortran translation of ranger 0.18.0 computational code.
! Upstream C++ core: Copyright (c) 2014-2018 Marvin N. Wright; MIT licensed.
! The upstream R package is GPL-3. See NOTICE.md for full provenance.
module ranger_maxstat
   use r_kinds, only : dp
   implicit none
   private

   public :: average_ranks, logrank_scores
   public :: maxstat_best, maxstat_pvalue_lau92, maxstat_pvalue_lau94
   public :: maxstat_pvalue_unadjusted, adjust_pvalues_bh

contains

   subroutine average_ranks(values, ranks)
      real(dp), intent(in) :: values(:)
      real(dp), intent(out) :: ranks(size(values))
      integer, allocatable :: order(:)
      integer :: i, j, n, reps, k
      real(dp) :: rank_value

      n = size(values)
      allocate(order(n))
      order = [(i, i = 1, n)]
      call sort_indices_real(values, order, decreasing=.false.)
      i = 1
      do while (i <= n)
         reps = 1
         do while (i + reps <= n)
            if (.not. real_equal(values(order(i + reps)), values(order(i)))) exit
            reps = reps + 1
         end do
         ! C++ ranger: (2*i0 + reps - 1)/2 + 1, with i0 = i - 1.
         rank_value = (2.0_dp * real(i - 1, dp) + real(reps, dp) - 1.0_dp) / 2.0_dp + 1.0_dp
         do j = 0, reps - 1
            k = order(i + j)
            ranks(k) = rank_value
         end do
         i = i + reps
      end do
   end subroutine average_ranks

   subroutine logrank_scores(time, status, scores)
      real(dp), intent(in) :: time(:)
      integer, intent(in) :: status(:)
      real(dp), intent(out) :: scores(size(time))
      integer, allocatable :: order(:)
      integer :: n, i, j, first_group
      real(dp) :: cumsum

      n = size(time)
      if (size(status) /= n) error stop 'logrank_scores: incompatible dimensions'
      allocate(order(n))
      order = [(i, i = 1, n)]
      call sort_indices_real(time, order, decreasing=.false.)
      scores = 0.0_dp
      cumsum = 0.0_dp
      first_group = 1
      do i = 1, n
         if (i < n) then
            if (real_equal(time(order(i)), time(order(i + 1)))) cycle
         end if
         do j = first_group, i
            cumsum = cumsum + real(status(order(j)), dp) / real(n - i + 1, dp)
         end do
         do j = first_group, i
            scores(order(j)) = real(status(order(j)), dp) - cumsum
         end do
         first_group = i + 1
      end do
   end subroutine logrank_scores

   subroutine maxstat_best(scores, x, minprop, maxprop, best_stat, best_split, num_left, found)
      real(dp), intent(in) :: scores(:), x(:), minprop, maxprop
      real(dp), intent(out) :: best_stat, best_split
      integer, allocatable, intent(out) :: num_left(:)
      logical, intent(out) :: found
      integer, allocatable :: order(:), temp_left(:)
      integer :: n, i, maxpos, minpos, nleft, ngroup, group_count
      real(dp) :: sum_all, mean_scores, sum_mean_diff, sum_scores, expected, variance, stat

      n = size(x)
      if (size(scores) /= n) error stop 'maxstat_best: incompatible dimensions'
      allocate(order(n), temp_left(n))
      order = [(i, i = 1, n)]
      call sort_indices_real(x, order, decreasing=.false.)

      sum_all = sum(scores)
      mean_scores = sum_all / real(n, dp)
      sum_mean_diff = sum((scores - mean_scores) ** 2)

      minpos = 1
      if (real(n, dp) * minprop > 1.0_dp) minpos = int(real(n, dp) * minprop - 1.0_dp) + 1
      maxpos = int(real(n, dp) * maxprop - 1.0_dp) + 1
      maxpos = min(n, max(1, maxpos))

      best_stat = -1.0_dp
      best_split = -1.0_dp
      sum_scores = 0.0_dp
      nleft = 0
      ngroup = 0
      group_count = 0
      do i = 1, n
         group_count = group_count + 1
         if (i < n) then
            if (real_equal(x(order(i)), x(order(i + 1)))) cycle
         end if
         ngroup = ngroup + 1
         if (ngroup == 1) then
            temp_left(ngroup) = group_count
         else
            temp_left(ngroup) = temp_left(ngroup - 1) + group_count
         end if
         group_count = 0
      end do
      allocate(num_left(ngroup))
      if (ngroup > 0) num_left = temp_left(1:ngroup)

      do i = 1, maxpos
         sum_scores = sum_scores + scores(order(i))
         nleft = nleft + 1
         if (i < minpos) cycle
         if (i < n) then
            if (real_equal(x(order(i)), x(order(i + 1)))) cycle
         end if
         if (real_equal(x(order(i)), x(order(n)))) exit
         variance = real(nleft * (n - nleft), dp) / real(n * (n - 1), dp) * sum_mean_diff
         if (variance <= 0.0_dp) cycle
         expected = real(nleft, dp) / real(n, dp) * sum_all
         stat = abs((sum_scores - expected) / sqrt(variance))
         if (stat > best_stat) then
            best_stat = stat
            if (i < n) then
               best_split = (x(order(i)) + x(order(i + 1))) / 2.0_dp
            else
               best_split = x(order(i))
            end if
         end if
      end do
      found = best_stat > -1.0_dp
   end subroutine maxstat_best

   pure real(dp) function maxstat_pvalue_lau92(b, minprop, maxprop) result(pvalue)
      real(dp), intent(in) :: b, minprop, maxprop
      real(dp) :: density, logprop

      if (b < 1.0_dp) then
         pvalue = 1.0_dp
         return
      end if
      logprop = log((maxprop * (1.0_dp - minprop)) / ((1.0_dp - maxprop) * minprop))
      density = exp(-0.5_dp * b * b) / sqrt(2.0_dp * acos(-1.0_dp))
      pvalue = 4.0_dp * density / b + density * (b - 1.0_dp / b) * logprop
      pvalue = max(pvalue, 0.0_dp)
   end function maxstat_pvalue_lau92

   pure real(dp) function maxstat_pvalue_lau94(b, minprop, maxprop, n, num_left) result(pvalue)
      real(dp), intent(in) :: b, minprop, maxprop
      integer, intent(in) :: n, num_left(:)
      real(dp) :: dterm, m1, m2, t, pi
      integer :: i

      ! minprop/maxprop are part of the upstream signature but Lau94 uses the
      ! realized admissible cutpoint counts directly.
      if (minprop < 0.0_dp .or. maxprop > 1.0_dp) error stop 'maxstat_pvalue_lau94: invalid proportions'
      pi = acos(-1.0_dp)
      dterm = 0.0_dp
      do i = 1, size(num_left) - 1
         m1 = real(num_left(i), dp)
         m2 = real(num_left(i + 1), dp)
         if (m1 <= 0.0_dp .or. m2 <= 0.0_dp .or. m1 >= real(n, dp)) cycle
         t = sqrt(max(0.0_dp, 1.0_dp - m1 * real(n - num_left(i + 1), dp) / &
            (real(n - num_left(i), dp) * m2)))
         dterm = dterm + exp(-0.5_dp * b * b) / pi * &
            (t - (b * b / 4.0_dp - 1.0_dp) * t ** 3 / 6.0_dp)
      end do
      pvalue = erfc(b / sqrt(2.0_dp)) + dterm
   end function maxstat_pvalue_lau94

   pure real(dp) function maxstat_pvalue_unadjusted(b) result(pvalue)
      real(dp), intent(in) :: b
      pvalue = erfc(b / sqrt(2.0_dp))
   end function maxstat_pvalue_unadjusted

   subroutine adjust_pvalues_bh(pvalue, adjusted)
      real(dp), intent(in) :: pvalue(:)
      real(dp), intent(out) :: adjusted(size(pvalue))
      integer, allocatable :: order(:)
      integer :: n, i, idx, last

      n = size(pvalue)
      if (n == 0) return
      allocate(order(n))
      order = [(i, i = 1, n)]
      call sort_indices_real(pvalue, order, decreasing=.true.)
      adjusted(order(1)) = pvalue(order(1))
      do i = 2, n
         idx = order(i)
         last = order(i - 1)
         adjusted(idx) = min(adjusted(last), real(n, dp) / real(n - i + 1, dp) * pvalue(idx))
      end do
   end subroutine adjust_pvalues_bh

   subroutine sort_indices_real(values, order, decreasing)
      real(dp), intent(in) :: values(:)
      integer, intent(inout) :: order(:)
      logical, intent(in) :: decreasing
      integer :: i, j, key

      do i = 2, size(order)
         key = order(i)
         j = i - 1
         if (decreasing) then
            do while (j >= 1)
               if (values(order(j)) >= values(key)) exit
               order(j + 1) = order(j)
               j = j - 1
            end do
         else
            do while (j >= 1)
               if (values(order(j)) <= values(key)) exit
               order(j + 1) = order(j)
               j = j - 1
            end do
         end if
         order(j + 1) = key
      end do
   end subroutine sort_indices_real

   pure logical function real_equal(a, b) result(equal)
      real(dp), intent(in) :: a, b
      equal = abs(a - b) <= 0.0_dp
   end function real_equal

end module ranger_maxstat
