! SPDX-License-Identifier: MIT
module r_quantiles
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   use r_ordering, only : r_order, r_sort_values_in_place
   implicit none
   private

   integer, parameter, public :: r_qrule_math = 1
   integer, parameter, public :: r_qrule_school = 2
   integer, parameter, public :: r_qrule_shahvaish = 3
   integer, parameter, public :: r_qrule_hf1 = 4
   integer, parameter, public :: r_qrule_hf2 = 5
   integer, parameter, public :: r_qrule_hf3 = 6
   integer, parameter, public :: r_qrule_hf4 = 7
   integer, parameter, public :: r_qrule_hf5 = 8
   integer, parameter, public :: r_qrule_hf6 = 9
   integer, parameter, public :: r_qrule_hf7 = 10
   integer, parameter, public :: r_qrule_hf8 = 11
   integer, parameter, public :: r_qrule_hf9 = 12

   public :: r_median, r_quantile_type7
   public :: r_weighted_quantile_ecdf, r_weighted_quantile_frequency_type7
   public :: r_weighted_quantile_isotone, r_weighted_quantile_linear_cdf
   public :: r_weighted_quantile_survey

contains

   pure real(dp) function r_median(x, na_rm) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: na_rm
      value = r_quantile_type7(x, 0.5_dp, na_rm)
   end function r_median

   pure real(dp) function r_quantile_type7(x, probability, na_rm) result(value)
      real(dp), intent(in) :: x(:), probability
      logical, intent(in), optional :: na_rm
      real(dp), allocatable :: values(:)
      real(dp) :: fraction, h
      integer :: j, n
      logical :: remove_na

      value = quiet_nan()
      if (.not. valid_probability(probability)) return
      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      if (.not. remove_na .and. any(ieee_is_nan(x))) return
      n = size(x)
      if (remove_na) n = count(.not. ieee_is_nan(x))
      if (n == 0) return
      allocate(values(n))
      if (remove_na) then
         values = pack(x, .not. ieee_is_nan(x))
      else
         values = x
      end if
      call sort_values(values)
      h = 1.0_dp + real(n - 1, dp)*probability
      j = int(floor(h))
      fraction = h - real(j, dp)
      if (j >= n) then
         value = values(n)
      else
         value = (1.0_dp - fraction)*values(j) + fraction*values(j + 1)
      end if
   end function r_quantile_type7

   pure real(dp) function r_weighted_quantile_ecdf(x, weights, probability, na_rm, finite_only) result(value)
      real(dp), intent(in) :: x(:), weights(:), probability
      logical, intent(in), optional :: na_rm, finite_only
      real(dp), allocatable :: xs(:), ws(:)
      real(dp) :: cumulative, target
      integer :: i
      logical :: ok

      call prepare_weighted(x, weights, probability, na_rm, finite_only, xs, ws, ok)
      value = quiet_nan()
      if (.not. ok) return
      target = probability*sum(ws)
      cumulative = 0.0_dp
      do i = 1, size(xs)
         cumulative = cumulative + ws(i)
         if (cumulative >= target) then
            value = xs(i)
            return
         end if
      end do
      value = xs(size(xs))
   end function r_weighted_quantile_ecdf

   pure real(dp) function r_weighted_quantile_linear_cdf(x, weights, probability, na_rm, finite_only) result(value)
      real(dp), intent(in) :: x(:), weights(:), probability
      logical, intent(in), optional :: na_rm, finite_only
      real(dp), allocatable :: cumulative(:), xs(:), ws(:)
      real(dp) :: fraction, total
      integer :: i
      logical :: ok

      call prepare_weighted(x, weights, probability, na_rm, finite_only, xs, ws, ok)
      value = quiet_nan()
      if (.not. ok) return
      allocate(cumulative(size(xs)))
      total = sum(ws)
      cumulative(1) = ws(1)/total
      do i = 2, size(xs)
         cumulative(i) = cumulative(i - 1) + ws(i)/total
      end do
      if (probability <= cumulative(1)) then
         value = xs(1)
         return
      end if
      do i = 2, size(xs)
         if (cumulative(i) > probability) then
            fraction = (probability - cumulative(i - 1))/(cumulative(i) - cumulative(i - 1))
            value = xs(i - 1) + fraction*(xs(i) - xs(i - 1))
            return
         end if
      end do
      value = xs(size(xs))
   end function r_weighted_quantile_linear_cdf

   pure real(dp) function r_weighted_quantile_frequency_type7(x, weights, probability, na_rm, finite_only) result(value)
      real(dp), intent(in) :: x(:), weights(:), probability
      logical, intent(in), optional :: na_rm, finite_only
      real(dp), allocatable :: xs(:), ws(:)
      real(dp) :: cumulative, fraction, target, total
      integer :: i
      logical :: ok

      call prepare_weighted(x, weights, probability, na_rm, finite_only, xs, ws, ok)
      value = quiet_nan()
      if (.not. ok) return
      total = sum(ws)
      target = 1.0_dp + (total - 1.0_dp)*probability
      cumulative = 0.0_dp
      do i = 1, size(xs)
         cumulative = cumulative + ws(i)
         if (cumulative >= target) then
            if (i == 1 .or. abs(cumulative - target) > 1.0e-12_dp .or. i == size(xs)) then
               value = xs(i)
            else
               fraction = target - floor(target)
               value = (1.0_dp - fraction)*xs(i) + fraction*xs(i + 1)
            end if
            return
         end if
      end do
      value = xs(size(xs))
   end function r_weighted_quantile_frequency_type7

   pure real(dp) function r_weighted_quantile_isotone(x, weights, probability, na_rm, finite_only) result(value)
      real(dp), intent(in) :: x(:), weights(:), probability
      logical, intent(in), optional :: na_rm, finite_only
      real(dp), allocatable :: xs(:), ws(:)
      real(dp) :: cumulative, target, tolerance, total
      integer :: i
      logical :: ok

      call prepare_weighted(x, weights, probability, na_rm, finite_only, xs, ws, ok)
      value = quiet_nan()
      if (.not. ok) return
      total = sum(ws)
      target = probability*total
      tolerance = 10.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(total))
      cumulative = 0.0_dp
      do i = 1, size(xs)
         cumulative = cumulative + ws(i)
         if (cumulative >= target - tolerance) then
            if (abs(cumulative - target) <= tolerance .and. i < size(xs)) then
               value = (ws(i)*xs(i) + ws(i + 1)*xs(i + 1))/(ws(i) + ws(i + 1))
            else
               value = xs(i)
            end if
            return
         end if
      end do
      value = xs(size(xs))
   end function r_weighted_quantile_isotone

   pure real(dp) function r_weighted_quantile_survey(x, weights, probability, rule, na_rm, finite_only) result(value)
      real(dp), intent(in) :: x(:), weights(:), probability
      integer, intent(in) :: rule
      logical, intent(in), optional :: na_rm, finite_only
      real(dp), allocatable :: cumulative(:), knots(:), unique_weights(:), unique_x(:), xs(:), ws(:)
      real(dp) :: fraction, lower_gap, mean_weight, total, upper_gap
      integer :: i, lower, next, n
      logical :: ok

      call prepare_weighted(x, weights, probability, na_rm, finite_only, xs, ws, ok)
      value = quiet_nan()
      if (.not. ok) return
      n = size(xs)
      if (n == 1) then
         value = xs(1)
         return
      end if
      total = sum(ws)
      select case (rule)
      case (r_qrule_math, r_qrule_hf1, r_qrule_school, r_qrule_hf2, r_qrule_hf4)
         cumulative = cumulative_sum(ws)
         lower = last_not_greater(cumulative, probability*total)
         next = min(n, lower + 1)
         lower_gap = probability - cumulative(lower)/total
         upper_gap = cumulative(next)/total - probability
         select case (rule)
         case (r_qrule_math, r_qrule_hf1)
            if (lower_gap <= 0.0_dp) then
               value = xs(lower)
            else
               value = xs(next)
            end if
         case (r_qrule_school, r_qrule_hf2)
            if (lower_gap <= 0.0_dp) then
               value = 0.5_dp*(xs(lower) + xs(next))
            else
               value = xs(next)
            end if
         case (r_qrule_hf4)
            if (abs(upper_gap + lower_gap) <= tiny(1.0_dp)) then
               value = xs(lower)
            else
               fraction = lower_gap/(upper_gap + lower_gap)
               value = xs(lower) + fraction*(xs(next) - xs(lower))
            end if
         end select
      case (r_qrule_hf3)
         call aggregate_ties(xs, ws, unique_x, unique_weights)
         cumulative = cumulative_sum(unique_weights)
         lower = last_not_greater(cumulative, probability*sum(unique_weights))
         next = min(size(unique_x), lower + 1)
         lower_gap = probability - cumulative(lower)/sum(unique_weights)
         if (lower_gap <= 0.0_dp .and. mod(lower, 2) == 0) then
            value = unique_x(lower)
         else
            value = unique_x(next)
         end if
      case (r_qrule_hf5)
         cumulative = cumulative_sum(ws)
         knots = (cumulative - 0.5_dp*ws)/total
         value = linear_interpolation(knots, xs, probability)
      case (r_qrule_hf6)
         cumulative = cumulative_sum(ws)
         knots = cumulative/(total + ws(n))
         value = linear_interpolation(knots, xs, probability)
      case (r_qrule_shahvaish)
         mean_weight = total/real(n, dp)
         cumulative = cumulative_sum(ws/mean_weight)
         knots = (cumulative + 0.5_dp - ws/(2.0_dp*mean_weight))/real(n + 1, dp)
         value = constant_interpolation(knots, xs, probability)
      case (r_qrule_hf7)
         cumulative = cumulative_sum(ws)
         allocate(knots(n))
         knots(1) = 0.0_dp
         knots(2:n) = cumulative(1:n - 1)/cumulative(n - 1)
         value = linear_interpolation(knots, xs, probability)
      case (r_qrule_hf8)
         cumulative = cumulative_sum(ws)
         allocate(knots(n))
         knots(1) = (2.0_dp/3.0_dp)*cumulative(1)/(total + ws(n)/3.0_dp)
         knots(2:n) = (cumulative(1:n - 1)/3.0_dp + 2.0_dp*cumulative(2:n)/3.0_dp)/(total + ws(n)/3.0_dp)
         value = linear_interpolation(knots, xs, probability)
      case (r_qrule_hf9)
         cumulative = cumulative_sum(ws)
         allocate(knots(n))
         knots(1) = (5.0_dp/8.0_dp)*cumulative(1)/(total + ws(n)/4.0_dp)
         knots(2:n) = (3.0_dp*cumulative(1:n - 1)/8.0_dp + 5.0_dp*cumulative(2:n)/8.0_dp)/(total + ws(n)/4.0_dp)
         value = linear_interpolation(knots, xs, probability)
      case default
         value = quiet_nan()
      end select
   end function r_weighted_quantile_survey

   pure subroutine prepare_weighted(x, weights, probability, na_rm, finite_only, xs, ws, ok)
      real(dp), intent(in) :: x(:), weights(:), probability
      logical, intent(in), optional :: na_rm, finite_only
      real(dp), allocatable, intent(out) :: xs(:), ws(:)
      logical, intent(out) :: ok
      integer, allocatable :: order(:)
      integer :: i, n
      logical :: remove_na, require_finite

      ok = .false.
      if (size(x) /= size(weights) .or. .not. valid_probability(probability)) return
      if (any(.not. ieee_is_finite(weights)) .or. any(weights < 0.0_dp)) return
      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      require_finite = .false.
      if (present(finite_only)) require_finite = finite_only
      n = 0
      do i = 1, size(x)
         if (weights(i) == 0.0_dp) cycle
         if (require_finite .and. .not. ieee_is_finite(x(i))) cycle
         if (remove_na .and. ieee_is_nan(x(i))) cycle
         if (ieee_is_nan(x(i))) return
         n = n + 1
      end do
      if (n == 0) return
      allocate(xs(n), ws(n))
      n = 0
      do i = 1, size(x)
         if (weights(i) == 0.0_dp) cycle
         if (require_finite .and. .not. ieee_is_finite(x(i))) cycle
         if (remove_na .and. ieee_is_nan(x(i))) cycle
         n = n + 1
         xs(n) = x(i)
         ws(n) = weights(i)
      end do
      call r_order(xs, order)
      xs = xs(order)
      ws = ws(order)
      ok = ieee_is_finite(sum(ws)) .and. sum(ws) > 0.0_dp
   end subroutine prepare_weighted

   pure logical function valid_probability(probability) result(valid)
      real(dp), intent(in) :: probability
      valid = ieee_is_finite(probability) .and. probability >= 0.0_dp .and. probability <= 1.0_dp
   end function valid_probability

   pure real(dp) function quiet_nan() result(value)
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

   pure function cumulative_sum(x) result(values)
      real(dp), intent(in) :: x(:)
      real(dp) :: values(size(x))
      integer :: i
      values(1) = x(1)
      do i = 2, size(x)
         values(i) = values(i - 1) + x(i)
      end do
   end function cumulative_sum

   pure integer function last_not_greater(values, target) result(position)
      real(dp), intent(in) :: values(:), target
      integer :: i
      position = 1
      do i = 1, size(values)
         if (values(i) <= target) position = i
      end do
   end function last_not_greater

   pure real(dp) function linear_interpolation(knots, values, probability) result(value)
      real(dp), intent(in) :: knots(:), values(:), probability
      integer :: i
      if (probability <= knots(1)) then
         value = values(1)
         return
      end if
      if (probability >= knots(size(knots))) then
         value = values(size(values))
         return
      end if
      do i = 1, size(knots) - 1
         if (probability >= knots(i) .and. probability <= knots(i + 1)) then
            if (abs(knots(i + 1) - knots(i)) <= epsilon(1.0_dp)*max(1.0_dp, abs(knots(i)))) then
               value = values(i + 1)
            else
               value = values(i) + (values(i + 1) - values(i))*(probability - knots(i))/(knots(i + 1) - knots(i))
            end if
            return
         end if
      end do
      value = values(size(values))
   end function linear_interpolation

   pure real(dp) function constant_interpolation(knots, values, probability) result(value)
      real(dp), intent(in) :: knots(:), values(:), probability
      integer :: i
      if (probability <= knots(1)) then
         value = values(1)
         return
      end if
      value = values(size(values))
      do i = 1, size(knots) - 1
         if (probability < knots(i + 1)) then
            value = values(i)
            return
         end if
      end do
   end function constant_interpolation

   pure subroutine aggregate_ties(x, weights, unique_x, unique_weights)
      real(dp), intent(in) :: x(:), weights(:)
      real(dp), allocatable, intent(out) :: unique_x(:), unique_weights(:)
      real(dp), allocatable :: work_weights(:), work_x(:)
      integer :: i, n
      allocate(work_x(size(x)), work_weights(size(x)))
      n = 0
      do i = 1, size(x)
         if (n == 0 .or. abs(x(i) - work_x(n)) > epsilon(1.0_dp)*max(1.0_dp, abs(x(i)), abs(work_x(n)))) then
            n = n + 1
            work_x(n) = x(i)
            work_weights(n) = weights(i)
         else
            work_weights(n) = work_weights(n) + weights(i)
         end if
      end do
      allocate(unique_x(n), unique_weights(n))
      unique_x = work_x(:n)
      unique_weights = work_weights(:n)
   end subroutine aggregate_ties

   pure subroutine sort_values(values)
      real(dp), intent(inout) :: values(:)
      if (size(values) < 2) return
      call r_sort_values_in_place(values)
   end subroutine sort_values

end module r_quantiles
