! SPDX-License-Identifier: GPL-2.0-or-later
! Translation of randtests 1.0.2 by Frederico Caeiro and Ayana Mateus.
module randtests_distributions
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, &
                                                                               ieee_quiet_nan, ieee_negative_inf
   use r_kinds, only: dp
   use r_special, only: r_log_choose
   use r_stability, only: r_log_sum_exp
   implicit none
   private
   integer, parameter, public :: max_runs_total = 100000
   integer, parameter, public :: max_exact_bartels_n = 9
   public :: druns, pruns, qruns, rruns, dbartelsrank, pbartelsrank
   public :: bartels_histogram
contains
   pure elemental real(dp) function druns(x, n1, n2, log_probability) result(value)
      real(dp), intent(in) :: x !! Number of runs; noninteger or infinite values have zero mass.
      integer, intent(in) :: n1 !! Positive first-class count, with n1+n2 <= max_runs_total.
      integer, intent(in) :: n2 !! Positive second-class count, with n1+n2 <= max_runs_total.
      logical, intent(in), optional :: log_probability !! Return natural log mass when true; default false.
      real(dp) :: a, b, log_mass
      integer :: k, r
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. valid_counts(n1, n2) .or. ieee_is_nan(x)) return
      log_mass = ieee_value(0.0_dp, ieee_negative_inf)
      if (ieee_is_finite(x)) then
         if (x >= 2.0_dp .and. x <= real(run_maximum(n1, n2), dp)) then
            r = int(x)
            if (x == real(r, dp)) then
               if (mod(r, 2) == 0) then
                  k = r / 2 - 1
                  log_mass = log(2.0_dp) + r_log_choose(n1 - 1, k) + r_log_choose(n2 - 1, k)
               else
                  k = (r - 1) / 2
                  a = r_log_choose(n1 - 1, k) + r_log_choose(n2 - 1, k - 1)
                  b = r_log_choose(n1 - 1, k - 1) + r_log_choose(n2 - 1, k)
                  log_mass = r_log_sum_exp([a, b])
               end if
               log_mass = log_mass - r_log_choose(n1 + n2, n1)
            end if
         end if
      end if
      value = exp(log_mass)
      if (present(log_probability)) then
         if (log_probability) value = log_mass
      end if
   end function druns

   pure elemental real(dp) function pruns(q, n1, n2, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: q !! CDF cutoff; lower tail is P(R <= floor(q)); infinities allowed.
      integer, intent(in) :: n1 !! Positive first-class count, with n1+n2 <= max_runs_total.
      integer, intent(in) :: n2 !! Positive second-class count, with n1+n2 <= max_runs_total.
      logical, intent(in), optional :: lower_tail !! Select P(R <= q) if true (default), otherwise P(R > q).
      logical, intent(in), optional :: log_probability !! Return natural log probability when true; default false.
      integer :: lo, hi, r, maximum
      real(dp) :: log_total, term
      logical :: lower
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. valid_counts(n1, n2) .or. ieee_is_nan(q)) return
      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      maximum = run_maximum(n1, n2)
      lo = 2
      hi = maximum
      if (q < 2.0_dp) then
         if (lower) hi = 1
      else if (q >= real(maximum, dp)) then
         if (.not. lower) lo = maximum + 1
      else if (lower) then
         hi = floor(q)
      else
         lo = floor(q) + 1
      end if
      log_total = ieee_value(0.0_dp, ieee_negative_inf)
      if (lo == 2 .and. hi == maximum) then
         log_total = 0.0_dp
      else
         do r = lo, hi
            term = druns(real(r, dp), n1, n2, .true.)
            if (r == lo) then
               log_total = term
            else
               log_total = r_log_sum_exp([log_total, term])
            end if
         end do
      end if
      log_total = min(0.0_dp, log_total)
      value = exp(log_total)
      if (present(log_probability)) then
         if (log_probability) value = log_total
      end if
   end function pruns

   pure elemental real(dp) function qruns(p, n1, n2, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: p !! Probability in [0,1], or its natural logarithm in [-Inf,0].
      integer, intent(in) :: n1 !! Positive first-class count, with n1+n2 <= max_runs_total.
      integer, intent(in) :: n2 !! Positive second-class count, with n1+n2 <= max_runs_total.
      logical, intent(in), optional :: lower_tail !! Interpret p as a lower-tail probability if true (default).
      logical, intent(in), optional :: log_probability !! Interpret p as a natural log probability when true.
      integer :: lo, hi, mid
      logical :: lower, log_p
      real(dp) :: log_target, tail, tolerance
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. valid_counts(n1, n2) .or. ieee_is_nan(p)) return
      lower = .true.
      log_p = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_probability)) log_p = log_probability
      if (log_p) then
         if (p > 0.0_dp) return
         log_target = p
      else
         if (p < 0.0_dp .or. p > 1.0_dp) return
         log_target = ieee_value(0.0_dp, ieee_negative_inf)
         if (p > 0.0_dp) log_target = log(p)
      end if
      lo = 2
      hi = run_maximum(n1, n2)
      if (log_target == 0.0_dp) then
         value = real(merge(hi, lo, lower), dp)
         return
      end if
      if (.not. ieee_is_finite(log_target)) then
         value = real(merge(lo, hi, lower), dp)
         return
      end if
      ! Discrete boundary comparisons tolerate a few ulps from exp/log round trips.
      tolerance = 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(log_target))
      ! Compare the requested tail directly, retaining tiny upper probabilities.
      do while (lo < hi)
         mid = lo + (hi - lo) / 2
         tail = pruns(real(mid, dp), n1, n2, lower, .true.)
         if (lower) then
            if (tail >= log_target - tolerance) then
               hi = mid
            else
               lo = mid + 1
            end if
         else
            if (tail <= log_target + tolerance) then
               hi = mid
            else
               lo = mid + 1
            end if
         end if
      end do
      value = real(lo, dp)
   end function qruns

   subroutine rruns(values, n1, n2, status)
      integer, intent(out) :: values(:) !! Caller-sized output of random run counts; set to -1 on invalid counts.
      integer, intent(in) :: n1 !! Positive first-class count, with n1+n2 <= max_runs_total.
      integer, intent(in) :: n2 !! Positive second-class count, with n1+n2 <= max_runs_total.
      integer, intent(out), optional :: status !! Zero on success, one if either count is invalid.
      real(dp) :: u
      integer :: i
      if (present(status)) status = 1
      values = -1
      if (.not. valid_counts(n1, n2)) return
      do i = 1, size(values)
         call random_number(u)
         values(i) = int(qruns(u, n1, n2))
      end do
      if (present(status)) status = 0
   end subroutine rruns

   pure elemental real(dp) function dbartelsrank(x, n, log_probability) result(value)
      real(dp), intent(in) :: x !! Sum of squared successive rank differences; nonintegers have zero mass.
      integer, intent(in) :: n !! Number of untied ranks, from 1 through max_exact_bartels_n (9).
      logical, intent(in), optional :: log_probability !! Return natural log mass when true; default false.
      integer, allocatable :: histogram(:)
      integer :: status, total
      real(dp) :: probability
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (ieee_is_nan(x)) return
      call bartels_histogram(n, histogram, status)
      if (status /= 0) return
      total = sum(histogram)
      probability = 0.0_dp
      if (ieee_is_finite(x)) then
         if (x >= 0.0_dp .and. x <= real(ubound(histogram, 1), dp)) then
            if (x == real(int(x), dp)) probability = real(histogram(int(x)), dp) / real(total, dp)
         end if
      end if
      value = probability_result(probability, log_probability)
   end function dbartelsrank

   pure elemental real(dp) function pbartelsrank(q, n, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: q !! Statistic cutoff; lower tail is P(NM <= q), including infinite cutoffs.
      integer, intent(in) :: n !! Number of untied ranks, from 1 through max_exact_bartels_n (9).
      logical, intent(in), optional :: lower_tail !! Select P(NM <= q) if true (default), otherwise P(NM > q).
      logical, intent(in), optional :: log_probability !! Return natural log probability when true; default false.
      integer, allocatable :: histogram(:)
      integer :: status, total, subtotal, i
      logical :: lower
      real(dp) :: probability
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (ieee_is_nan(q)) return
      call bartels_histogram(n, histogram, status)
      if (status /= 0) return
      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      total = sum(histogram)
      subtotal = 0
      do i = 0, ubound(histogram, 1)
         if (lower .eqv. (real(i, dp) <= q)) subtotal = subtotal + histogram(i)
      end do
      probability = real(subtotal, dp) / real(total, dp)
      value = probability_result(probability, log_probability)
   end function pbartelsrank

   pure subroutine bartels_histogram(n, histogram, status)
      integer, intent(in) :: n !! Number of distinct ranks; exact enumeration supports 1 <= n <= 9.
      integer, allocatable, intent(out) :: histogram(:) !! Counts indexed from 0 by squared-difference statistic.
      integer, intent(out) :: status !! Zero on success; one for n outside the supported exact range.
      logical, allocatable :: used(:)
      status = 1
      if (n < 1 .or. n > max_exact_bartels_n) then
         allocate (histogram(0:-1))
         return
      end if
      allocate (histogram(0:n * (n - 1) ** 2), used(n))
      histogram = 0
      used = .false.
      call enumerate_bartels(n, 1, 0, 0, used, histogram)
      status = 0
   end subroutine bartels_histogram

   pure recursive subroutine enumerate_bartels(n, depth, previous, score, used, histogram)
      integer, intent(in) :: n !! Total number of distinct integer ranks.
      integer, intent(in) :: depth !! Current permutation position, starting at one.
      integer, intent(in) :: previous !! Rank chosen at the preceding position; zero at the root.
      integer, intent(in) :: score !! Squared-difference sum accumulated before the current position.
      logical, intent(inout) :: used(:) !! Length-n mask recording ranks already in the current permutation.
      integer, intent(inout) :: histogram(0:) !! Exact permutation counts indexed by completed score.
      integer :: i, next_score
      if (depth > n) then
         histogram(score) = histogram(score) + 1
         return
      end if
      do i = 1, n
         if (used(i)) cycle
         next_score = score
         if (depth > 1) next_score = score + (i - previous) ** 2
         used(i) = .true.
         call enumerate_bartels(n, depth + 1, i, next_score, used, histogram)
         used(i) = .false.
      end do
   end subroutine enumerate_bartels

   pure elemental logical function valid_counts(n1, n2)
      integer, intent(in) :: n1 !! First-class count checked for positivity and resource bounds.
      integer, intent(in) :: n2 !! Second-class count checked before integer addition can overflow.
      valid_counts = .false.
      if (n1 < 1 .or. n2 < 1) return
      if (n1 > max_runs_total .or. n2 > max_runs_total) return
      valid_counts = n1 + n2 <= max_runs_total
   end function valid_counts

   pure elemental integer function run_maximum(n1, n2)
      integer, intent(in) :: n1 !! Valid positive first-class count.
      integer, intent(in) :: n2 !! Valid positive second-class count.
      run_maximum = 2 * min(n1, n2)
      if (n1 /= n2) run_maximum = run_maximum + 1
   end function run_maximum

   pure elemental real(dp) function probability_result(p, return_log) result(value)
      real(dp), intent(in) :: p !! Probability in [0,1] produced by exact enumeration.
      logical, intent(in), optional :: return_log !! Return log probability if true, with log(0) represented as -Inf.
      value = p
      if (.not. present(return_log)) return
      if (.not. return_log) return
      value = ieee_value(0.0_dp, ieee_negative_inf)
      if (p > 0.0_dp) value = log(p)
   end function probability_result
end module randtests_distributions
