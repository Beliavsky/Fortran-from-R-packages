! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from A.I. McLeod's Kendall 2.2.2. See NOTICE.md for provenance.
module kendall_module
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use r_ordering, only : r_sort_values_in_place
   use r_distributions, only : r_pnorm
   implicit none
   private

   integer, parameter, public :: kendall_ok = 0
   integer, parameter, public :: kendall_size_mismatch = 1
   integer, parameter, public :: kendall_too_few = 2
   integer, parameter, public :: kendall_nonfinite = 3
   integer, parameter, public :: kendall_degenerate = 4
   integer, parameter, public :: kendall_invalid_method = 5
   integer, parameter, public :: kendall_exact_ties = 6
   integer, parameter, public :: kendall_exact_limit = 7
   integer, parameter, public :: kendall_invalid_period = 8
   integer, parameter, public :: kendall_asymptotic = 1
   integer, parameter, public :: kendall_exact = 2
   integer, parameter, public :: kendall_exact_max_n = 100

   type, public :: kendall_result
      !! Five numerical fields corresponding to the R Kendall result, plus diagnostics.
      real(dp) :: tau = 0.0_dp  !! Tau-b, or NaN if undefined or the request is invalid.
      real(dp) :: sl = 0.0_dp   !! Two-sided p-value, or NaN when not computed.
      real(dp) :: s = 0.0_dp    !! Concordant minus discordant pair count.
      real(dp) :: d = 0.0_dp    !! Tie-corrected denominator of tau-b.
      real(dp) :: vars = 0.0_dp !! Null variance of the score with ties accounted for.
      integer :: n_used = 0     !! Number of complete observations after NaN omission.
      integer :: status = kendall_ok !! Zero on success, otherwise a public status constant.
      logical :: exact = .false. !! True only when the optional no-ties exact method was used.
   end type kendall_result

   public :: dp, kendall, mann_kendall, seasonal_mann_kendall

contains

   pure subroutine kendall(x, y, result, method)
      !! Compute Kendall tau-b and its two-sided significance level.
      !! Default: the continuity-corrected normal approximation used by the upstream runtime.
      real(dp), intent(in) :: x(:) !! First variable; pairs with either value NaN are omitted.
      real(dp), intent(in) :: y(:) !! Second variable, same length as x; retained infinities are invalid.
      type(kendall_result), intent(out) :: result !! Statistics, complete-pair count, and diagnostic status.
      integer, intent(in), optional :: method !! kendall_asymptotic (default) or kendall_exact (no ties, n <= 100).
      real(dp), allocatable :: a(:), b(:)
      real(dp) :: tx, ax, bx, ty, ay, by, rn, cn, pairs, adjusted
      integer :: n, i, j, selected

      call initialize_result(result)
      selected = kendall_asymptotic
      if (present(method)) selected = method
      if (selected /= kendall_asymptotic .and. selected /= kendall_exact) then
         result%status = kendall_invalid_method
         return
      end if
      if (size(x) /= size(y)) then
         result%status = kendall_size_mismatch
         return
      end if
      n = count(.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
      result%n_used = n
      if (n < 3) then
         result%status = kendall_too_few
         return
      end if
      allocate(a(n), b(n))
      j = 0
      do i = 1, size(x)
         if (ieee_is_nan(x(i)) .or. ieee_is_nan(y(i))) cycle
         if (.not. ieee_is_finite(x(i)) .or. .not. ieee_is_finite(y(i))) then
            result%status = kendall_nonfinite
            return
         end if
         j = j + 1
         a(j) = x(i)
         b(j) = y(i)
      end do

      result%s = 0.0_dp
      do i = 1, n - 1
         do j = i + 1, n
            result%s = result%s + pair_score(a(i), b(i), a(j), b(j))
         end do
      end do
      call tie_moments(a, tx, ax, bx)
      call tie_moments(b, ty, ay, by)
      rn = real(n, dp)
      cn = rn * (rn - 1.0_dp)
      pairs = 0.5_dp * cn
      result%d = sqrt(pairs - 0.5_dp * tx) * sqrt(pairs - 0.5_dp * ty)
      if (tx == 0.0_dp .and. ty == 0.0_dp) result%d = pairs
      if (result%d == 0.0_dp) then
         result%vars = 0.0_dp
         result%status = kendall_degenerate
         return
      end if
      result%tau = max(-1.0_dp, min(1.0_dp, result%s / result%d))
      result%vars = (cn * (2.0_dp * rn + 5.0_dp) - ax - ay) / 18.0_dp &
         + bx * by / (9.0_dp * cn * (rn - 2.0_dp)) + tx * ty / (2.0_dp * cn)

      if (selected == kendall_exact) then
         if (tx > 0.0_dp .or. ty > 0.0_dp) then
            result%status = kendall_exact_ties
            return
         end if
         if (n > kendall_exact_max_n) then
            result%status = kendall_exact_limit
            return
         end if
         result%sl = exact_probability(n, nint(abs(result%s)))
         result%exact = .true.
      else
         adjusted = max(0.0_dp, abs(result%s) - 1.0_dp)
         result%sl = 2.0_dp * r_pnorm(adjusted / sqrt(result%vars), lower_tail = .false.)
      end if
   end subroutine kendall

   pure subroutine mann_kendall(x, result, method)
      !! Test monotonic trend against observation positions in the original input.
      real(dp), intent(in) :: x(:) !! Ordered observations; NaNs are omitted without reordering other values.
      type(kendall_result), intent(out) :: result !! Kendall statistics and number of complete observations.
      integer, intent(in), optional :: method !! kendall_asymptotic (default) or no-ties kendall_exact.
      real(dp), allocatable :: time(:)
      integer :: i

      time = [(real(i, dp), i = 1, size(x))]
      call kendall(time, x, result, method)
   end subroutine mann_kendall

   pure subroutine seasonal_mann_kendall(x, period, result)
      !! Aggregate within-season scores, denominators and variances as in SeasonalMannKendall.
      !! The final normal significance level has no continuity correction, matching upstream.
      real(dp), intent(in) :: x(:) !! Consecutive seasonal observations; NaNs preserve their time slots.
      integer, intent(in) :: period !! Positive seasons per cycle; every season needs at least three non-NaNs.
      type(kendall_result), intent(out) :: result !! Aggregated statistics; constant seasons contribute zeros.
      type(kendall_result) :: one_season
      real(dp) :: score, denom, variance
      integer :: season, total

      call initialize_result(result)
      result%n_used = count(.not. ieee_is_nan(x))
      if (period < 1) then
         result%status = kendall_invalid_period
         return
      end if
      if (period > size(x) / 3) then
         result%status = kendall_too_few
         return
      end if
      score = 0.0_dp
      denom = 0.0_dp
      variance = 0.0_dp
      total = 0
      do season = 1, period
         call mann_kendall(x(season::period), one_season)
         if (one_season%status /= kendall_ok .and. one_season%status /= kendall_degenerate) then
            result%status = one_season%status
            return
         end if
         total = total + one_season%n_used
         score = score + one_season%s
         denom = denom + one_season%d
         variance = variance + one_season%vars
      end do
      result%n_used = total
      result%s = score
      result%d = denom
      result%vars = variance
      if (denom == 0.0_dp .or. variance == 0.0_dp) then
         result%status = kendall_degenerate
         return
      end if
      result%tau = score / denom
      result%sl = 2.0_dp * r_pnorm(abs(score) / sqrt(variance), lower_tail = .false.)
   end subroutine seasonal_mann_kendall

   pure subroutine initialize_result(result)
      type(kendall_result), intent(out) :: result !! Result reset to NaNs, zero count, success, and non-exact method.

      result%tau = ieee_value(0.0_dp, ieee_quiet_nan)
      result%sl = result%tau
      result%s = result%tau
      result%d = result%tau
      result%vars = result%tau
      result%n_used = 0
      result%status = kendall_ok
      result%exact = .false.
   end subroutine initialize_result

   pure elemental real(dp) function pair_score(x1, y1, x2, y2) result(score)
      real(dp), intent(in) :: x1 !! Finite first-variable observation in the first pair position.
      real(dp), intent(in) :: y1 !! Finite second-variable observation in the first pair position.
      real(dp), intent(in) :: x2 !! Finite first-variable observation in the second pair position.
      real(dp), intent(in) :: y2 !! Finite second-variable observation in the second pair position.

      score = 0.0_dp
      if (x1 == x2 .or. y1 == y2) return
      if ((x1 < x2 .and. y1 < y2) .or. (x1 > x2 .and. y1 > y2)) then
         score = 1.0_dp
      else
         score = -1.0_dp
      end if
   end function pair_score

   pure subroutine tie_moments(values, pairs, variance, triples)
      real(dp), intent(inout) :: values(:) !! Finite complete observations, sorted in-place for tie grouping.
      real(dp), intent(out) :: pairs !! Sum over tie groups of t times (t - 1), including ordered pairs.
      real(dp), intent(out) :: variance !! Sum over tie groups of t times (t - 1) times (2t + 5).
      real(dp), intent(out) :: triples !! Sum over tie groups of t times (t - 1) times (t - 2).
      real(dp) :: t
      integer :: first, last, n

      call r_sort_values_in_place(values)
      pairs = 0.0_dp
      variance = 0.0_dp
      triples = 0.0_dp
      n = size(values)
      first = 1
      do while (first <= n)
         last = first + 1
         do while (last <= n)
            if (values(last) /= values(first)) exit
            last = last + 1
         end do
         t = real(last - first, dp)
         pairs = pairs + t * (t - 1.0_dp)
         variance = variance + t * (t - 1.0_dp) * (2.0_dp * t + 5.0_dp)
         triples = triples + t * (t - 1.0_dp) * (t - 2.0_dp)
         first = last
      end do
   end subroutine tie_moments

   pure real(dp) function exact_probability(n, abs_score) result(probability)
      !! Independent no-ties extension, not a translation of the unreachable upstream AS 71 branch.
      !! Inversion-count probabilities follow the polynomial product of (1 + z + ... + z**(m - 1))/m.
      integer, intent(in) :: n !! Validated no-ties sample size, between 3 and kendall_exact_max_n.
      integer, intent(in) :: abs_score !! Absolute concordance score, with valid permutation parity.
      real(dp), allocatable :: previous(:), current(:)
      integer :: cutoff, m, k, top, old_top, lower, upper

      cutoff = (n * (n - 1) / 2 - abs_score) / 2
      allocate(previous(0:cutoff), current(0:cutoff))
      previous = 0.0_dp
      previous(0) = 1.0_dp
      old_top = 0
      do m = 2, n
         top = min(cutoff, m * (m - 1) / 2)
         current = 0.0_dp
         do k = 0, top
            lower = max(0, k - m + 1)
            upper = min(k, old_top)
            current(k) = sum(previous(lower:upper)) / real(m, dp)
         end do
         previous = current
         old_top = top
      end do
      probability = min(1.0_dp, 2.0_dp * sum(previous))
   end function exact_probability

end module kendall_module
