! SPDX-License-Identifier: Artistic-2.0 AND GPL-2.0-or-later
! Numerical translation of matrixStats 1.5.0 ordering and weighted summaries.
! Original matrixStats implementations: Henrik Bengtsson (2007-2015),
! Robert Gentleman, Hector Corrada Bravo, Harris Jaffee, Peter Langfelder,
! Brian Montgomery, and Ola Hossjer. See NOTICE and upstream provenance.
module matrixstats_order
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf
   use r_kinds, only: dp
   use r_ordering, only: r_order, r_sort_values_in_place
   implicit none
   private
   public :: ms_quantile, ms_median, ms_iqr, ms_mad, ms_ranks, ms_order_stat
   public :: ms_weighted_mean, ms_weighted_median, ms_weighted_var, ms_weighted_sd, ms_weighted_mad
contains
   pure function ms_quantile(x, probability, na_rm, type) result(value)
      !! Hyndman-Fan sample quantiles, with matrixStats missing-value semantics.
      real(dp), intent(in) :: x(:) !! Observations; NaNs denote missing observations.
      real(dp), intent(in) :: probability !! Quantile probability in [0,1], allowing 100 machine epsilons of rounding.
      logical, intent(in), optional :: na_rm !! Drop NaNs before sorting; defaults to false.
      integer, intent(in), optional :: type !! Hyndman-Fan estimator number 1 through 9; defaults to 7.
      real(dp) :: value, p, h, fraction, fuzz
      real(dp), allocatable :: values(:)
      integer :: estimator, n, j, lower, upper
      logical :: ok
      value = missing_value()
      if (.not. ieee_is_finite(probability)) return
      if (probability < -100.0_dp * epsilon(1.0_dp)) return
      if (probability > 1.0_dp + 100.0_dp * epsilon(1.0_dp)) return
      estimator = 7
      if (present(type)) estimator = type
      if (estimator < 1 .or. estimator > 9) return
      call nonmissing_values(x, na_rm, values, ok)
      if (.not. ok) return
      n = size(values)
      if (n == 0) return
      call r_sort_values_in_place(values)
      p = max(0.0_dp, min(1.0_dp, probability))
      select case (estimator)
      case (1, 2, 4)
         h = real(n, dp) * p
      case (3)
         h = real(n, dp) * p - 0.5_dp
      case (5)
         h = real(n, dp) * p + 0.5_dp
      case (6)
         h = real(n + 1, dp) * p
      case (7)
         h = 1.0_dp + real(n - 1, dp) * p
      case (8)
         h = (real(n, dp) + 1.0_dp / 3.0_dp) * p + 1.0_dp / 3.0_dp
      case (9)
         h = (real(n, dp) + 0.25_dp) * p + 0.375_dp
      end select
      fuzz = 4.0_dp * epsilon(1.0_dp)
      if (estimator == 7) fuzz = 0.0_dp
      j = floor(h * (1.0_dp + fuzz))
      fraction = h - real(j, dp)
      select case (estimator)
      case (1)
         fraction = merge(1.0_dp, 0.0_dp, fraction > 0.0_dp)
      case (2)
         fraction = merge(1.0_dp, 0.5_dp, fraction > 0.0_dp)
      case (3)
         fraction = merge(1.0_dp, 0.0_dp, fraction /= 0.0_dp .or. mod(j, 2) /= 0)
      case (4:9)
         if (abs(fraction) < fuzz) fraction = 0.0_dp
      end select
      lower = max(1, min(n, j))
      upper = max(1, min(n, j + 1))
      if (fraction <= 0.0_dp) then
         value = values(lower)
      else if (fraction >= 1.0_dp) then
         value = values(upper)
      else if (values(lower) == values(upper)) then
         value = values(lower)
      else
         value = (1.0_dp - fraction) * values(lower) + fraction * values(upper)
      end if
   end function ms_quantile

   pure function ms_median(x, na_rm) result(value)
      !! Median using matrixStats' arithmetic average for an even sample size.
      real(dp), intent(in) :: x(:) !! Observations; NaNs denote missing observations.
      logical, intent(in), optional :: na_rm !! Drop NaNs before sorting; defaults to false.
      real(dp) :: value
      real(dp), allocatable :: values(:)
      integer :: n
      logical :: ok
      value = missing_value()
      call nonmissing_values(x, na_rm, values, ok)
      if (.not. ok) return
      n = size(values)
      if (n == 0) return
      call r_sort_values_in_place(values)
      if (mod(n, 2) == 1) then
         value = values(n / 2 + 1)
      else
         value = (values(n / 2) + values(n / 2 + 1)) / 2.0_dp
      end if
   end function ms_median

   pure function ms_iqr(x, na_rm, type) result(value)
      !! Scalar iqr(), including its zero result for a single nonmissing value.
      real(dp), intent(in) :: x(:) !! Observations; NaNs denote missing observations.
      logical, intent(in), optional :: na_rm !! Drop NaNs before computing quartiles; defaults to false.
      integer, intent(in), optional :: type !! Hyndman-Fan estimator number 1 through 9; defaults to 7.
      real(dp) :: value
      real(dp), allocatable :: values(:)
      logical :: ok
      value = missing_value()
      call nonmissing_values(x, na_rm, values, ok)
      if (.not. ok) return
      if (size(values) == 0) return
      if (size(values) == 1) then
         value = 0.0_dp
      else
         value = ms_quantile(values, 0.75_dp, type=type) - ms_quantile(values, 0.25_dp, type=type)
      end if
   end function ms_iqr

   pure function ms_mad(x, na_rm, constant, center) result(value)
      !! Median absolute deviation, as used by rowMads and colMads.
      real(dp), intent(in) :: x(:) !! Observations; NaNs denote missing observations.
      logical, intent(in), optional :: na_rm !! Drop original NaNs before computing the MAD; defaults to false.
      real(dp), intent(in), optional :: constant !! Multiplicative scale factor; defaults to 1.4826.
      real(dp), intent(in), optional :: center !! Fixed center; absent means use the sample median.
      real(dp) :: value, location, scale
      real(dp), allocatable :: values(:)
      logical :: ok
      value = missing_value()
      call nonmissing_values(x, na_rm, values, ok)
      if (.not. ok) return
      if (size(values) == 0) return
      scale = 1.4826_dp
      if (present(constant)) scale = constant
      if (present(center)) then
         location = center
      else
         if (size(values) == 1) then
            value = 0.0_dp
            return
         end if
         location = ms_median(values)
      end if
      if (present(center)) then
         value = scale * ms_median(abs(values - location), na_rm=na_rm)
      else
         value = scale * ms_median(abs(values - location))
      end if
   end function ms_mad

   pure function ms_order_stat(x, which) result(value)
      !! One-based order statistic; invalid indices and missing observations return NaN.
      real(dp), intent(in) :: x(:) !! Observations to order; any NaN makes the result NaN.
      integer, intent(in) :: which !! One-based index in the sorted observations, between 1 and size(x).
      real(dp) :: value
      real(dp), allocatable :: values(:)
      value = missing_value()
      if (which < 1 .or. which > size(x)) return
      if (any(ieee_is_nan(x))) return
      values = x
      call r_sort_values_in_place(values)
      value = values(which)
   end function ms_order_stat

   pure function ms_ranks(x, ties_method, na_last) result(ranks)
      !! Deterministic rank methods, including dense ranks; random ranks are not supported.
      real(dp), intent(in) :: x(:) !! Observations ranked in ascending order; infinities are valid values.
      character(len=*), intent(in), optional :: ties_method !! max (default), min, average, first, last, or dense.
      character(len=*), intent(in), optional :: na_last !! keep (default), first, or last; keep returns NaN for missing entries.
      real(dp) :: ranks(size(x))
      real(dp), allocatable :: values(:)
      integer, allocatable :: indices(:), order(:)
      character(len=:), allocatable :: method, missing
      integer :: i, j, k, n, offset, dense, nmissing
      ranks = missing_value()
      method = 'max'
      if (present(ties_method)) method = ties_method
      select case (method)
      case ('max', 'min', 'average', 'first', 'last', 'dense')
      case default
         return
      end select
      missing = 'keep'
      if (present(na_last)) missing = na_last
      select case (missing)
      case ('keep', 'first', 'last')
      case default
         return
      end select
      indices = pack([(i, i = 1, size(x))], .not. ieee_is_nan(x))
      values = x(indices)
      n = size(values)
      nmissing = size(x) - n
      offset = 0
      if (missing == 'first') offset = nmissing
      call r_order(values, order)
      i = 1
      dense = 0
      do while (i <= n)
         j = i
         do while (j < n)
            if (values(order(j + 1)) /= values(order(i))) exit
            j = j + 1
         end do
         dense = dense + 1
         do k = i, j
            select case (method)
            case ('min')
               ranks(indices(order(k))) = real(i + offset, dp)
            case ('max')
               ranks(indices(order(k))) = real(j + offset, dp)
            case ('average')
               ranks(indices(order(k))) = real(i + j, dp) / 2.0_dp + real(offset, dp)
            case ('first')
               ranks(indices(order(k))) = real(k + offset, dp)
            case ('last')
               ranks(indices(order(k))) = real(i + j - k + offset, dp)
            case ('dense')
               ranks(indices(order(k))) = real(dense + offset, dp)
            end select
         end do
         i = j + 1
      end do
      if (missing /= 'keep') then
         offset = 0
         if (missing == 'last') offset = n
         do i = 1, size(x)
            if (.not. ieee_is_nan(x(i))) cycle
            offset = offset + 1
            ranks(i) = real(offset, dp)
         end do
      end if
   end function ms_ranks

   pure function ms_weighted_mean(x, w, na_rm, refine) result(value)
      !! weightedMean: zero weights are omitted, but negative weights contribute as upstream does.
      real(dp), intent(in) :: x(:) !! Observations; NaNs are missing and zero-weight observations are ignored.
      real(dp), intent(in), optional :: w(:) !! Weights matching x; absent means unit weights. Negative weights are retained.
      logical, intent(in), optional :: na_rm !! Omit pairs with missing x; missing weights otherwise propagate. Default false.
      logical, intent(in), optional :: refine !! Apply a residual correction to a finite mean; defaults to false.
      real(dp) :: value, weight, total, wtotal, residual, correction, wcorrection, term
      integer :: i
      logical :: remove, improve
      value = missing_value()
      if (present(w)) then
         if (size(w) /= size(x)) return
      end if
      remove = .false.
      if (present(na_rm)) remove = na_rm
      improve = .false.
      if (present(refine)) improve = refine
      total = 0.0_dp
      wtotal = 0.0_dp
      correction = 0.0_dp
      wcorrection = 0.0_dp
      do i = 1, size(x)
         weight = 1.0_dp
         if (present(w)) weight = w(i)
         if (weight == 0.0_dp) cycle
         if (ieee_is_nan(x(i))) then
            if (remove) cycle
            return
         end if
         if (ieee_is_nan(weight)) return
         term = weight * x(i)
         call compensated_add(total, correction, term)
         call compensated_add(wtotal, wcorrection, weight)
      end do
      total = total + correction
      wtotal = wtotal + wcorrection
      if (.not. ieee_is_nan(wtotal) .and. .not. ieee_is_finite(wtotal)) return
      if (.not. ieee_is_finite(total)) then
         value = total
         return
      end if
      if (wtotal == 0.0_dp) then
         if (ieee_is_nan(total) .or. total == 0.0_dp) return
         value = sign(ieee_value(0.0_dp, ieee_positive_inf), total * sign(1.0_dp, wtotal))
         return
      end if
      value = total / wtotal
      if (.not. improve .or. .not. ieee_is_finite(value)) return
      residual = 0.0_dp
      correction = 0.0_dp
      do i = 1, size(x)
         weight = 1.0_dp
         if (present(w)) weight = w(i)
         if (weight == 0.0_dp) cycle
         if (ieee_is_nan(x(i)) .and. remove) cycle
         call compensated_add(residual, correction, weight * (x(i) - value))
      end do
      value = value + (residual + correction) / wtotal
   end function ms_weighted_mean

   pure function ms_weighted_median(x, w, na_rm, interpolate, ties) result(value)
      !! matrixStats midpoint-CDF interpolation or a classical weighted median.
      real(dp), intent(in) :: x(:) !! Observations; zero and negative weights normally exclude observations.
      real(dp), intent(in), optional :: w(:) !! Weights matching x; absent means unit weights. Infinite weights dominate.
      logical, intent(in), optional :: na_rm !! Omit missing observations and missing weights as upstream C does; default false.
      logical, intent(in), optional :: interpolate !! Interpolate midpoint cumulative weights; default true unless ties is given.
      character(len=*), intent(in), optional :: ties !! Classical ties: weighted (default), min, max, or mean.
      real(dp) :: value, total, cumulative, previous, midpoint, weight, low, high
      real(dp), allocatable :: values(:), weights(:), active_weights(:)
      integer, allocatable :: order(:)
      character(len=:), allocatable :: method
      integer :: i, j, n
      logical :: remove, linear, infinite_weights
      value = missing_value()
      if (present(w)) then
         if (size(w) /= size(x)) return
      end if
      remove = .false.
      if (present(na_rm)) remove = na_rm
      linear = .not. present(ties)
      if (present(interpolate)) linear = interpolate
      method = 'weighted'
      if (present(ties)) method = ties
      select case (method)
      case ('weighted', 'min', 'max', 'mean')
      case default
         return
      end select
      allocate(active_weights(size(x)))
      active_weights = 0.0_dp
      infinite_weights = .false.
      ! Preserve the upstream C filter order, including the early +Inf branch.
      do i = 1, size(x)
         weight = 1.0_dp
         if (present(w)) weight = w(i)
         if (ieee_is_nan(weight)) then
            if (.not. remove) return
         else if (weight <= 0.0_dp) then
            cycle
         else if (.not. ieee_is_finite(weight)) then
            infinite_weights = .true.
            active_weights = 0.0_dp
            do j = 1, size(x)
               weight = w(j)
               if (ieee_is_nan(weight)) then
                  if (.not. remove) return
               else if (.not. ieee_is_finite(weight)) then
                  if (ieee_is_nan(x(j))) then
                     if (.not. remove) return
                  else
                     active_weights(j) = 1.0_dp
                  end if
               end if
            end do
            exit
         else if (ieee_is_nan(x(i))) then
            if (.not. remove) return
         else
            active_weights(i) = weight
         end if
      end do
      values = pack(x, active_weights > 0.0_dp)
      weights = pack(active_weights, active_weights > 0.0_dp)
      n = size(values)
      if (n == 0) return
      if (n == 1) then
         value = values(1)
         return
      end if
      if (infinite_weights) then
         ! Upstream ignores ties/interpolate when infinite weights occur.
         ! Return the actual median, repairing its odd-size partial-sort off-by-one bug.
         value = ms_median(values)
         return
      end if
      call weighted_median_order(values, order)
      weights = weights(order)
      total = sum(weights)
      if (.not. ieee_is_finite(total)) return
      cumulative = 0.0_dp
      previous = 0.0_dp
      do i = 1, n
         weight = weights(i) / total
         cumulative = cumulative + weight
         if (linear) then
            midpoint = cumulative - weight / 2.0_dp
            if (midpoint >= 0.5_dp) then
               if (i == 1) then
                  value = values(i)
               else
                  value = values(i) + (0.5_dp - midpoint) / (midpoint - previous) * (values(i) - values(i - 1))
               end if
               return
            end if
            previous = midpoint
         else
            if (cumulative > 0.5_dp) then
               if (i == 1) then
                  value = values(i)
                  return
               end if
               low = previous
               high = 1.0_dp - low
               if (high > 0.5_dp) then
                  value = values(i)
               else
                  select case (method)
                  case ('min')
                     value = values(i - 1)
                  case ('max')
                     value = values(i)
                  case ('mean')
                     value = (values(i - 1) + values(i)) / 2.0_dp
                  case ('weighted')
                     value = low * values(i - 1) + high * values(i)
                  end select
               end if
               return
            end if
            previous = cumulative
         end if
      end do
      value = values(n)
   end function ms_weighted_median

   pure function ms_weighted_var(x, w, na_rm, center) result(value)
      !! Frequency-weight variance with denominator sum(w)-1, rather than reliability-weight correction.
      real(dp), intent(in) :: x(:) !! Observations; zero or negative weights remove the corresponding observations.
      real(dp), intent(in), optional :: w(:) !! Nonnegative frequency weights matching x; absent means unit weights.
      logical, intent(in), optional :: na_rm !! Remove pairs with missing x, but retain other missing weights; default false.
      real(dp), intent(in), optional :: center !! Fixed center; absent means the weighted arithmetic mean.
      real(dp) :: value, location, total, numerator
      real(dp), allocatable :: values(:), weights(:)
      logical :: ok
      value = missing_value()
      call prepare_positive_weights(x, w, na_rm, values, weights, ok)
      if (.not. ok) return
      if (size(values) <= 1) return
      total = sum(weights)
      if (present(center)) then
         location = center
      else
         location = sum(weights * values) / total
      end if
      numerator = sum(weights * (values - location)**2)
      if (total == 1.0_dp) then
         if (ieee_is_nan(numerator) .or. numerator == 0.0_dp) return
         value = ieee_value(0.0_dp, ieee_positive_inf)
      else
         value = numerator / (total - 1.0_dp)
      end if
   end function ms_weighted_var

   pure function ms_weighted_sd(x, w, na_rm, center) result(value)
      !! Square root of the frequency-weight variance; negative variance gives NaN.
      real(dp), intent(in) :: x(:) !! Observations; zero or negative weights remove the corresponding observations.
      real(dp), intent(in), optional :: w(:) !! Frequency weights matching x; absent means unit weights.
      logical, intent(in), optional :: na_rm !! Remove pairs with missing x, but retain other missing weights; default false.
      real(dp), intent(in), optional :: center !! Fixed center; absent means the weighted arithmetic mean.
      real(dp) :: value
      value = ms_weighted_var(x, w, na_rm, center)
      if (value < 0.0_dp) then
         value = missing_value()
      else
         value = sqrt(value)
      end if
   end function ms_weighted_sd

   pure function ms_weighted_mad(x, w, na_rm, constant, center) result(value)
      !! Weighted median absolute deviation, using the upstream interpolated weighted median.
      real(dp), intent(in) :: x(:) !! Observations; zero and negative weights remove observations when weights are supplied.
      real(dp), intent(in), optional :: w(:) !! Weights matching x; absent uses the ordinary median absolute deviation.
      logical, intent(in), optional :: na_rm !! Remove pairs with missing x, but retain other missing weights; default false.
      real(dp), intent(in), optional :: constant !! Multiplicative normal-consistency factor; defaults to 1.4826.
      real(dp), intent(in), optional :: center !! Fixed center; absent means the corresponding weighted median.
      real(dp) :: value, location, scale
      real(dp), allocatable :: values(:), weights(:)
      logical :: ok
      value = missing_value()
      scale = 1.4826_dp
      if (present(constant)) scale = constant
      if (.not. present(w)) then
         call nonmissing_values(x, na_rm, values, ok)
         if (.not. ok) return
         if (present(center)) then
            location = center
         else
            location = ordinary_median(values)
         end if
         value = scale * ordinary_median(abs(values - location))
         return
      end if
      call prepare_positive_weights(x, w, na_rm, values, weights, ok)
      if (.not. ok) return
      if (size(values) == 0) return
      if (size(values) == 1) then
         value = 0.0_dp
         return
      end if
      if (present(center)) then
         location = center
      else
         location = ms_weighted_median(values, weights)
      end if
      value = scale * ms_weighted_median(abs(values - location), weights)
   end function ms_weighted_mad

   pure function ordinary_median(x) result(value)
      !! stats::median fallback for unweighted weightedMad, avoiding overflow in an even median.
      real(dp), intent(in) :: x(:) !! Already filtered observations or residuals; remaining NaNs propagate.
      real(dp) :: value, left, right
      real(dp), allocatable :: values(:)
      integer :: n
      value = missing_value()
      n = size(x)
      if (n == 0) return
      if (any(ieee_is_nan(x))) return
      values = x
      call r_sort_values_in_place(values)
      if (mod(n, 2) == 1) then
         value = values(n / 2 + 1)
      else
         left = values(n / 2)
         right = values(n / 2 + 1)
         if (ieee_is_finite(left) .and. ieee_is_finite(right)) then
            if (abs(left) > huge(1.0_dp) / 2.0_dp .or. abs(right) > huge(1.0_dp) / 2.0_dp) then
               value = left / 2.0_dp + right / 2.0_dp
               return
            end if
         end if
         value = (left + right) / 2.0_dp
      end if
   end function ordinary_median

   pure subroutine nonmissing_values(x, na_rm, values, ok)
      real(dp), intent(in) :: x(:) !! Raw observations, with IEEE NaNs treated as missing.
      logical, intent(in), optional :: na_rm !! Whether to discard NaNs; defaults to false.
      real(dp), allocatable, intent(out) :: values(:) !! Retained observations in original order, or an empty array on failure.
      logical, intent(out) :: ok !! True unless an unremoved missing observation was found.
      logical :: remove
      remove = .false.
      if (present(na_rm)) remove = na_rm
      ok = remove .or. .not. any(ieee_is_nan(x))
      if (ok) then
         values = pack(x, .not. ieee_is_nan(x))
      else
         allocate(values(0))
      end if
   end subroutine nonmissing_values

   pure subroutine prepare_positive_weights(x, w, na_rm, values, weights, ok)
      real(dp), intent(in) :: x(:) !! Raw observations before filtering nonpositive weights and missing observations.
      real(dp), intent(in), optional :: w(:) !! Raw weights matching x; absent means unit weights.
      logical, intent(in), optional :: na_rm !! Whether to remove pairs with missing observations; defaults to false.
      real(dp), allocatable, intent(out) :: values(:) !! Retained observations in original order.
      real(dp), allocatable, intent(out) :: weights(:) !! Positive weights, or unit weights if a retained weight is infinite.
      logical, intent(out) :: ok !! False for incompatible lengths or remaining missing observations or weights.
      real(dp), allocatable :: all_weights(:)
      logical, allocatable :: keep(:)
      logical :: remove
      ok = .false.
      allocate(values(0), weights(0))
      if (present(w)) then
         if (size(w) /= size(x)) return
         all_weights = w
      else
         allocate(all_weights(size(x)))
         all_weights = 1.0_dp
      end if
      keep = ieee_is_nan(all_weights) .or. all_weights > 0.0_dp
      remove = .false.
      if (present(na_rm)) remove = na_rm
      if (remove) keep = keep .and. .not. ieee_is_nan(x)
      values = pack(x, keep)
      weights = pack(all_weights, keep)
      if (any(ieee_is_nan(values)) .or. any(ieee_is_nan(weights))) return
      if (any(.not. ieee_is_finite(weights))) then
         values = pack(values, .not. ieee_is_finite(weights))
         weights = pack(weights, .not. ieee_is_finite(weights))
         weights = 1.0_dp
      end if
      ok = .true.
   end subroutine prepare_positive_weights

   pure subroutine weighted_median_order(values, order)
      !! Exact R_qsort_I permutation, needed because weighted interpolation depends on tie order.
      !! Shared r_order is stable and remains in use for deterministic ranking.
      ! R : A Computer Language for Statistical Data Analysis
      ! Copyright (C) 2002-2012 The R Core Team.
      ! Based on CACM algorithm #347 by R. C. Singleton (1969).
      !
      ! This program is free software; you can redistribute it and/or modify
      ! it under the terms of the GNU General Public License as published by
      ! the Free Software Foundation; either version 2 of the License, or
      ! (at your option) any later version.
      !
      ! This program is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! along with this program; if not, a copy is available at
      ! https://www.R-project.org/Licenses/
      !
      ! Translated from src/main/qsort-body.c, R-4-5-branch commit
      ! b4322f3207e8d99d9dc6e6fa6b80ad23c4dc2f5f. The control-flow labels
      ! retain the exact unstable tie permutation of R's public C entry point.
      real(dp), intent(inout) :: values(:) !! Nonmissing observations, sorted ascending on return.
      integer, allocatable, intent(out) :: order(:) !! Original one-based positions after R's unstable quicksort.
      integer :: lower_stack(bit_size(0) + 1), upper_stack(bit_size(0) + 1)
      integer :: i, j, k, l, m, pivot_index, pivot_order, saved_order
      real(dp) :: pivot, saved, fraction
      order = [(i, i = 1, size(values))]
      if (size(values) <= 1) return
      i = 1
      j = size(values)
      m = 1
      fraction = 0.375_dp
10    continue
      if (i >= j) go to 80
      if (fraction < 0.5898437_dp) then
         fraction = fraction + 0.0390625_dp
      else
         fraction = fraction - 0.21875_dp
      end if
20    continue
      k = i
      pivot_index = i + int(real(j - i, dp) * fraction)
      pivot_order = order(pivot_index)
      pivot = values(pivot_index)
      if (values(i) > pivot) then
         order(pivot_index) = order(i)
         order(i) = pivot_order
         pivot_order = order(pivot_index)
         values(pivot_index) = values(i)
         values(i) = pivot
         pivot = values(pivot_index)
      end if
      l = j
      if (values(j) < pivot) then
         order(pivot_index) = order(j)
         order(j) = pivot_order
         pivot_order = order(pivot_index)
         values(pivot_index) = values(j)
         values(j) = pivot
         pivot = values(pivot_index)
         if (values(i) > pivot) then
            order(pivot_index) = order(i)
            order(i) = pivot_order
            pivot_order = order(pivot_index)
            values(pivot_index) = values(i)
            values(i) = pivot
            pivot = values(pivot_index)
         end if
      end if
      do
         do
            l = l - 1
            if (values(l) <= pivot) exit
         end do
         saved_order = order(l)
         saved = values(l)
         do
            k = k + 1
            if (values(k) >= pivot) exit
         end do
         if (k > l) exit
         order(l) = order(k)
         order(k) = saved_order
         values(l) = values(k)
         values(k) = saved
      end do
      m = m + 1
      if (l - i <= j - k) then
         lower_stack(m) = k
         upper_stack(m) = j
         j = l
      else
         lower_stack(m) = i
         upper_stack(m) = l
         i = k
      end if
      go to 90
80    continue
      if (m == 1) return
      i = lower_stack(m)
      j = upper_stack(m)
      m = m - 1
90    continue
      if (j - i > 10) go to 20
      if (i == 1) go to 10
      i = i - 1
100   continue
      do
         i = i + 1
         if (i == j) go to 80
         pivot_order = order(i + 1)
         pivot = values(i + 1)
         if (values(i) > pivot) exit
      end do
      k = i
      do
         order(k + 1) = order(k)
         values(k + 1) = values(k)
         k = k - 1
         if (pivot >= values(k)) exit
      end do
      order(k + 1) = pivot_order
      values(k + 1) = pivot
      go to 100
   end subroutine weighted_median_order

   pure subroutine compensated_add(total, correction, term)
      real(dp), intent(inout) :: total !! Running uncorrected sum, including any encountered infinity.
      real(dp), intent(inout) :: correction !! Neumaier residual accumulated for finite additions.
      real(dp), intent(in) :: term !! Next summand; NaNs and infinities follow IEEE arithmetic.
      real(dp) :: updated
      updated = total + term
      if (ieee_is_finite(updated) .and. ieee_is_finite(total) .and. ieee_is_finite(term)) then
         if (abs(total) >= abs(term)) then
            correction = correction + ((total - updated) + term)
         else
            correction = correction + ((term - updated) + total)
         end if
      else
         correction = 0.0_dp
      end if
      total = updated
   end subroutine compensated_add

   pure function missing_value() result(value)
      real(dp) :: value
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function missing_value
end module matrixstats_order
