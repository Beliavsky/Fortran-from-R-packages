! SPDX-License-Identifier: GPL-2.0-only
! The matrixTests KS interface follows R/kolmogorov.R by Karolis Koncevicius.
! The integer-boundary probability recursion below is an original implementation.
module matrixtests_ks
   !! Two-sample Kolmogorov-Smirnov statistics and exact or limiting p-values.
   !! Exact probabilities condition on the pooled observations, including ties.
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use r_ordering, only : r_order
   implicit none
   private

   integer, parameter, public :: ks_success = 0
   integer, parameter, public :: ks_invalid_sample = 1
   integer, parameter, public :: ks_invalid_alternative = 2
   integer, parameter, public :: ks_resource_limit = 3
   integer(int64), parameter, public :: ks_max_exact_cells = 25000000_int64
   public :: ks_test

contains

   pure subroutine ks_test(x, y, statistic, pvalue, status, alternative, exact, used_exact)
      !! Test equality of two distributions, using the ECDF of x minus that of y.
      !! Positive ECDF excursions are tested by "greater", negative by "less".
      !! The exact method is a direct upper-tail random-label recursion and never
      !! falls back to an approximation. Its n*m work limit is 25,000,000 cells.
      real(dp), intent(in) :: x(:) !! Nonempty first sample, unsorted, allowing infinities but no NaNs.
      real(dp), intent(in) :: y(:) !! Nonempty second sample, unsorted, allowing infinities but no NaNs.
      real(dp), intent(out) :: statistic !! Maximum requested ECDF deviation in [0,1], or NaN on invalid input.
      real(dp), intent(out) :: pvalue !! Inclusive upper-tail probability, or NaN when status is nonzero.
      integer, intent(out) :: status !! 0 success, 1 invalid sample, 2 invalid alternative, or 3 resource limit.
      character(*), intent(in), optional :: alternative !! "two.sided" (default), "less", "greater", or a unique prefix.
      logical, intent(in), optional :: exact !! Force exact or asymptotic probabilities, defaulting to exact when n*m < 10000.
      logical, intent(out), optional :: used_exact !! True only when an exact probability was successfully evaluated.
      real(dp), allocatable :: pooled(:)
      integer, allocatable :: order(:)
      logical, allocatable :: boundary(:)
      integer :: n, m, total, k, nx, ny, mode, allocation_status
      integer(int64) :: product, deviation, threshold
      real(dp) :: effective_size
      logical :: do_exact, end_group

      statistic = ieee_value(0.0_dp, ieee_quiet_nan)
      pvalue = statistic
      status = ks_invalid_sample
      if (present(used_exact)) used_exact = .false.
      n = size(x)
      m = size(y)
      if (n == 0 .or. m == 0) return
      if (any(ieee_is_nan(x)) .or. any(ieee_is_nan(y))) return
      mode = 0
      if (present(alternative)) then
         if (len_trim(alternative) == 0) then
            status = ks_invalid_alternative
            return
         else if (index("two.sided", trim(alternative)) == 1) then
            mode = 0
         else if (index("greater", trim(alternative)) == 1) then
            mode = 1
         else if (index("less", trim(alternative)) == 1) then
            mode = -1
         else
            status = ks_invalid_alternative
            return
         end if
      end if

      product = int(n, int64) * int(m, int64)
      do_exact = product < 10000_int64
      if (present(exact)) do_exact = exact
      status = ks_resource_limit
      if (n > huge(total) - m) return
      if (do_exact .and. product > ks_max_exact_cells) return
      total = n + m
      allocate(pooled(total), boundary(0:total), stat = allocation_status)
      if (allocation_status /= 0) return
      pooled(:n) = x
      pooled(n + 1:) = y
      call r_order(pooled, order)
      boundary = .false.
      nx = 0
      ny = 0
      threshold = 0_int64
      do k = 1, total
         if (order(k) <= n) then
            nx = nx + 1
         else
            ny = ny + 1
         end if
         end_group = k == total
         if (k < total) end_group = pooled(order(k)) /= pooled(order(k + 1))
         if (.not. end_group) cycle
         boundary(k) = .true.
         deviation = int(nx, int64) * int(m, int64) - int(ny, int64) * int(n, int64)
         threshold = max(threshold, selected_deviation(deviation, mode))
      end do
      statistic = real(threshold, dp) / real(product, dp)
      if (do_exact) then
         call exact_upper_tail(n, m, boundary, threshold, mode, pvalue, status)
         if (present(used_exact)) used_exact = status == ks_success
      else
         effective_size = real(n, dp) * (real(m, dp) / real(total, dp))
         if (mode == 0) then
            pvalue = kolmogorov_upper_tail(sqrt(effective_size) * statistic)
         else
            pvalue = exp(-2.0_dp * effective_size * statistic**2)
         end if
         status = ks_success
      end if
   end subroutine ks_test

   pure elemental function selected_deviation(deviation, mode) result(value)
      integer(int64), intent(in) :: deviation !! Integer-scaled ECDF difference, x count times m minus y count times n.
      integer, intent(in) :: mode !! 0 for absolute deviation, 1 for positive deviation, or -1 for negative deviation.
      integer(int64) :: value

      select case (mode)
      case (0)
         value = abs(deviation)
      case (1)
         value = deviation
      case default
         value = -deviation
      end select
   end function selected_deviation

   pure subroutine exact_upper_tail(n, m, boundary, threshold, mode, probability, status)
      !! At state (i,j), average the probabilities at the two possible next labels.
      !! The remaining-label probabilities generate the uniform conditional law
      !! over all choose(n+m,n) assignments. A path hits the rejection region only
      !! at the end of a tied group. An absorbing hit has probability one, avoiding
      !! subtracting a nearly-unit CDF from one for extreme upper tails.
      integer, intent(in) :: n !! Positive size of the first sample.
      integer, intent(in) :: m !! Positive size of the second sample.
      logical, intent(in) :: boundary(0:) !! True at pooled ranks ending a distinct-value group, length n+m+1.
      integer(int64), intent(in) :: threshold !! Observed nonnegative statistic multiplied by n*m exactly.
      integer, intent(in) :: mode !! 0 two-sided, 1 greater, or -1 less, referring to the ECDF difference x-y.
      real(dp), intent(out) :: probability !! Inclusive conditional upper-tail probability without CDF subtraction.
      integer, intent(out) :: status !! 0 on success or 3 when the working-array allocation fails.
      real(dp), allocatable :: next(:)
      real(dp) :: remaining, value
      integer(int64) :: deviation
      integer :: i, j, allocation_status

      status = ks_success
      if (threshold == 0_int64) then
         probability = 1.0_dp
         return
      end if
      allocate(next(0:m), stat = allocation_status)
      if (allocation_status /= 0) then
         status = ks_resource_limit
         probability = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      next = 0.0_dp
      do i = n, 0, -1
         do j = m, 0, -1
            if (boundary(i + j)) then
               deviation = int(i, int64) * int(m, int64) - int(j, int64) * int(n, int64)
               if (selected_deviation(deviation, mode) >= threshold) then
                  next(j) = 1.0_dp
                  cycle
               end if
            end if
            if (i == n .and. j == m) then
               next(j) = 0.0_dp
               cycle
            end if
            remaining = real(n - i + m - j, dp)
            value = 0.0_dp
            if (i < n) value = real(n - i, dp) / remaining * next(j)
            if (j < m) value = value + real(m - j, dp) / remaining * next(j + 1)
            next(j) = value
         end do
      end do
      probability = min(1.0_dp, max(0.0_dp, next(0)))
   end subroutine exact_upper_tail

   pure function kolmogorov_upper_tail(z) result(probability)
      !! Evaluate the complementary Kolmogorov limit distribution. The Jacobi
      !! theta transformation is used for z<1 and the alternating survival series
      !! for z>=1, so neither slow small-z convergence nor tail cancellation occurs.
      real(dp), intent(in) :: z !! Nonnegative ECDF maximum multiplied by sqrt(n*m/(n+m)).
      real(dp) :: probability
      real(dp) :: total, term, exponent, sign, pi, odd
      integer :: k

      if (z <= 0.0_dp) then
         probability = 1.0_dp
         return
      end if
      total = 0.0_dp
      if (z < 1.0_dp) then
         pi = acos(-1.0_dp)
         exponent = -(pi / z)**2 / 8.0_dp
         do k = 1, 1000
            odd = real(2 * k - 1, dp)
            term = exp(exponent * odd**2)
            total = total + term
            if (term <= epsilon(total) * total) exit
         end do
         probability = 1.0_dp - sqrt(2.0_dp * pi) / z * total
      else
         exponent = -2.0_dp * z**2
         sign = 1.0_dp
         do k = 1, 1000
            term = exp(exponent * real(k, dp)**2)
            total = total + sign * term
            if (term <= epsilon(total) * abs(total)) exit
            sign = -sign
         end do
         probability = 2.0_dp * total
      end if
      probability = min(1.0_dp, max(0.0_dp, probability))
   end function kolmogorov_upper_tail

end module matrixtests_ks
