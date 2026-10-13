! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation of outliers 0.15 by Lukasz Komsta.
module outliers_distributions
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use r_distributions, only : r_qt, r_pt, r_qf, r_pf
   use outliers_tables
   implicit none
   private
   public :: qtable, qdixon, pdixon, qgrubbs, pgrubbs, qcochran, pcochran
   public :: dixon_type
   interface qtable
      module procedure qtable_scalar, qtable_vector
   end interface
contains
   pure real(dp) function qtable_scalar(p, probs, quants) result(value)
      real(dp), intent(in) :: p !! Evaluation coordinate; finite extrapolation is permitted.
      real(dp), intent(in) :: probs(:) !! Finite tabulated coordinates; order is arbitrary and local duplicates follow R.
      real(dp), intent(in) :: quants(:) !! Corresponding finite table values, same length as probs.
      real(dp) :: x(size(probs)), y(size(probs)), dist(size(probs)), a, b, term
      integer :: i, j, k, m, first, second, left
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      m = size(probs)
      if (m /= size(quants) .or. m < 2 .or. .not. ieee_is_finite(p)) return
      if (any(.not. ieee_is_finite(probs)) .or. any(.not. ieee_is_finite(quants))) return
      x = probs
      y = quants
      do i = 2, m
         a = x(i)
         b = y(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= a) exit
            x(j + 1) = x(j)
            y(j + 1) = y(j)
            j = j - 1
         end do
         x(j + 1) = a
         y(j + 1) = b
      end do
      if (p <= x(1) .or. p >= x(m)) then
         first = 1
         if (p >= x(m)) first = m - 1
         if (x(first + 1) == x(first)) then
            value = 0.5_dp * (y(first) + y(first + 1))
         else
            value = y(first) + (p - x(first)) * (y(first + 1) - y(first)) / (x(first + 1) - x(first))
         end if
         return
      end if
      if (m < 4) return
      dist = abs(p - x)
      first = minloc(dist, dim = 1)
      dist(first) = huge(1.0_dp)
      second = minloc(dist, dim = 1)
      left = max(2, min(min(first, second), m - 2)) - 1
      if (any(x(left + 1:left + 3) <= x(left:left + 2))) return
      value = 0.0_dp
      do j = left, left + 3
         term = y(j)
         do k = left, left + 3
            if (k == j) cycle
            term = term * (p - x(k)) / (x(j) - x(k))
         end do
         value = value + term
      end do
   end function qtable_scalar

   pure function qtable_vector(p, probs, quants) result(values)
      real(dp), intent(in) :: p(:) !! Finite evaluation coordinates, returned in the same order.
      real(dp), intent(in) :: probs(:) !! Finite tabulated coordinates; local duplicates follow R.
      real(dp), intent(in) :: quants(:) !! Corresponding finite table values, same length as probs.
      real(dp) :: values(size(p))
      integer :: i
      do i = 1, size(p)
         values(i) = qtable_scalar(p(i), probs, quants)
      end do
   end function qtable_vector

   pure elemental integer function dixon_type(n, variant) result(kind)
      integer, intent(in) :: n !! Number of nonmissing observations, at most 30.
      integer, intent(in) :: variant !! Dixon ratio code 0, 10, 11, 12, 20, 21, or 22; 0 selects by n.
      integer :: minimum
      kind = variant
      if (kind == 0) then
         if (n <= 7) then
            kind = 10
         else if (n <= 10) then
            kind = 11
         else if (n <= 13) then
            kind = 21
         else
            kind = 22
         end if
      end if
      select case (kind)
      case (10)
         minimum = 3
      case (11, 20)
         minimum = 4
      case (12, 21)
         minimum = 5
      case (22)
         minimum = 6
      case default
         kind = -1
         return
      end select
      if (n < minimum .or. n > 30) kind = -1
   end function dixon_type

   pure elemental real(dp) function qdixon(p, n, type, rev) result(value)
      real(dp), intent(in) :: p !! Probability, or Dixon ratio when rev is true; extrapolation is clamped to [0,1].
      integer, intent(in) :: n !! Sample size in the upstream type-specific range ending at 30.
      integer, intent(in), optional :: type !! Dixon code 0, 10 (default), 11, 12, 20, 21, or 22.
      logical, intent(in), optional :: rev !! Reverse the table interpolation to obtain a p-value; default false.
      integer :: kind, index
      logical :: reverse
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      kind = 10
      if (present(type)) kind = type
      kind = dixon_type(n, kind)
      if (kind < 0 .or. .not. ieee_is_finite(p)) return
      reverse = .false.
      if (present(rev)) reverse = rev
      index = 3 * (kind / 10 - 1) + mod(kind, 10) + 1
      if (reverse) then
         value = qtable(p, dixon_tables(:, n, index), dixon_probabilities)
      else
         value = qtable(p, dixon_probabilities, dixon_tables(:, n, index))
      end if
      if (.not. ieee_is_nan(value)) value = min(1.0_dp, max(0.0_dp, value))
   end function qdixon

   pure elemental real(dp) function pdixon(q, n, type) result(value)
      real(dp), intent(in) :: q !! Dixon ratio; finite values outside the table use endpoint extrapolation.
      integer, intent(in) :: n !! Sample size in the upstream type-specific range ending at 30.
      integer, intent(in), optional :: type !! Dixon ratio code; default 10.
      value = qdixon(q, n, type, rev = .true.)
   end function pdixon

   pure elemental real(dp) function qgrubbs(p, n, type, rev) result(value)
      real(dp), intent(in) :: p !! Probability, or test statistic when rev is true.
      integer, intent(in) :: n !! Sample size at least 3; type 20 requires 4 through 30.
      integer, intent(in), optional :: type !! Grubbs variant 10 (default), 11, or 20.
      logical, intent(in), optional :: rev !! Return the corresponding distribution function; default false.
      real(dp) :: rn, t, denominator, lower, upper, bound, probability, coordinate, tolerance
      integer :: kind
      logical :: reverse
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      kind = 10
      if (present(type)) kind = type
      reverse = .false.
      if (present(rev)) reverse = rev
      if (n < 3 .or. .not. ieee_is_finite(p)) return
      rn = real(n, dp)
      select case (kind)
      case (10, 11)
         if (.not. reverse) then
            if (p < 0.0_dp .or. p >= 1.0_dp) return
            probability = (1.0_dp - p) / rn
            if (kind == 11) probability = probability / (rn - 1.0_dp)
            t = r_qt(probability, rn - 2.0_dp)
            if (kind == 10) then
               value = (rn - 1.0_dp) / sqrt(rn) / sqrt(1.0_dp + (rn - 2.0_dp) / t ** 2)
            else
               value = sqrt(2.0_dp * (rn - 1.0_dp) / (1.0_dp + (rn - 2.0_dp) / t ** 2))
            end if
         else if (kind == 10) then
            denominator = (rn - 1.0_dp) ** 2 - p ** 2 * rn
            if (denominator <= 0.0_dp) then
               value = 1.0_dp
            else
               t = abs(p) * sqrt(rn * (rn - 2.0_dp) / denominator)
               value = 1.0_dp - min(1.0_dp, rn * r_pt(t, rn - 2.0_dp, lower_tail = .false.))
            end if
         else
            ! Preserve upstream clipping and the uniroot lower-bracket gap.
            t = r_qt((1.0_dp - 0.9999_dp) / (rn * (rn - 1.0_dp)), rn - 2.0_dp)
            upper = sqrt(2.0_dp * (rn - 1.0_dp) / (1.0_dp + (rn - 2.0_dp) / t ** 2))
            t = r_qt((1.0_dp - 2.0e-16_dp) / (rn * (rn - 1.0_dp)), rn - 2.0_dp)
            lower = sqrt(2.0_dp * (rn - 1.0_dp) / (1.0_dp + (rn - 2.0_dp) / t ** 2))
            t = r_qt((1.0_dp - 0.001_dp) / (rn * (rn - 1.0_dp)), rn - 2.0_dp)
            bound = sqrt(2.0_dp * (rn - 1.0_dp) / (1.0_dp + (rn - 2.0_dp) / t ** 2))
            ! Treat roundoff-sized boundary differences as equality across t backends.
            coordinate = p
            tolerance = 128.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(p))
            if (abs(p - lower) <= tolerance) coordinate = lower
            if (abs(p - bound) <= tolerance) coordinate = bound
            if (abs(p - upper) <= tolerance) coordinate = upper
            if (coordinate > upper) then
               value = 1.0_dp
            else if (coordinate < lower) then
               value = 0.0_dp
            else if (coordinate >= bound) then
               denominator = 2.0_dp * (rn - 1.0_dp) - coordinate ** 2
               t = coordinate * sqrt((rn - 2.0_dp) / denominator)
               value = 1.0_dp - rn * (rn - 1.0_dp) * r_pt(t, rn - 2.0_dp, lower_tail = .false.)
            end if
         end if
      case (20)
         if (n < 4 .or. n > 30) return
         if (reverse) then
            value = qtable(p, grubbs_table(:, n), grubbs_probabilities)
         else
            value = qtable(p, grubbs_probabilities, grubbs_table(:, n))
         end if
         if (.not. ieee_is_nan(value)) value = min(1.0_dp, max(0.0_dp, value))
      end select
   end function qgrubbs

   pure elemental real(dp) function pgrubbs(q, n, type) result(value)
      real(dp), intent(in) :: q !! Grubbs G statistic (types 10/11), or U ratio (type 20).
      integer, intent(in) :: n !! Sample size at least 3; type 20 requires 4 through 30.
      integer, intent(in), optional :: type !! Grubbs variant 10 (default), 11, or 20.
      value = qgrubbs(q, n, type, rev = .true.)
   end function pgrubbs

   pure elemental real(dp) function pcochran(q, n, k) result(value)
      real(dp), intent(in) :: q !! Variance ratio C, in [0,1].
      real(dp), intent(in) :: n !! Mean number of observations per group, greater than 1.
      integer, intent(in) :: k !! Number of groups, at least 2.
      real(dp) :: rk, f
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. ieee_is_finite(q) .or. .not. ieee_is_finite(n)) return
      if (q < 0.0_dp .or. q > 1.0_dp .or. n <= 1.0_dp .or. k < 2) return
      value = 0.0_dp
      if (q == 0.0_dp) return
      rk = real(k, dp)
      f = (1.0_dp / q - 1.0_dp) / (rk - 1.0_dp)
      value = min(1.0_dp, max(0.0_dp, 1.0_dp - rk * r_pf(f, (n - 1.0_dp) * (rk - 1.0_dp), n - 1.0_dp)))
   end function pcochran

   pure elemental real(dp) function qcochran(p, n, k) result(value)
      real(dp), intent(in) :: p !! Cumulative probability in [0,1].
      real(dp), intent(in) :: n !! Mean number of observations per group, greater than 1.
      integer, intent(in) :: k !! Number of groups, at least 2.
      real(dp) :: rk, f
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. ieee_is_finite(p) .or. .not. ieee_is_finite(n)) return
      if (p < 0.0_dp .or. p > 1.0_dp .or. n <= 1.0_dp .or. k < 2) return
      rk = real(k, dp)
      f = r_qf((1.0_dp - p) / rk, (n - 1.0_dp) * (rk - 1.0_dp), n - 1.0_dp)
      value = 1.0_dp / (1.0_dp + (rk - 1.0_dp) * f)
   end function qcochran
end module outliers_distributions
