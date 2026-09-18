! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the R package roll by Jason Foster; see PROVENANCE.md and NOTICE.md.
module roll_univariate
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: roll_na_logical = -huge(0)
   public :: roll_any_vec, roll_all_vec, roll_sum_vec, roll_prod_vec, roll_mean_vec
   public :: roll_min_vec, roll_max_vec, roll_idxmin_vec, roll_idxmax_vec
   public :: roll_median_vec, roll_quantile_vec, roll_var_vec, roll_sd_vec, roll_scale_vec
   public :: roll_pair_vec

contains

   pure function na_real() result(value)
      real(dp) :: value
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function na_real

   pure integer function effective_min_obs(width, min_obs) result(value)
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      integer, intent(in), optional :: min_obs !! Minimum finite observations; defaults to `width`.
      value = width
      if (present(min_obs)) value = min_obs
   end function effective_min_obs

   pure real(dp) function window_weight(weights, count) result(value)
      real(dp), intent(in), optional :: weights(:) !! Observation weights; the last weight aligns with the current row.
      integer, intent(in) :: count !! Zero-based lag within the current window.
      if (present(weights)) then
         value = weights(size(weights) - count)
      else
         value = 1.0_dp
      end if
   end function window_weight

   pure subroutine roll_any_vec(x, width, result, min_obs, complete_obs, na_restore, online)
      integer, intent(in) :: x(:) !! Logical series encoded as zero/one with `roll_na_logical` for missing values.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      integer, allocatable, intent(out) :: result(:) !! Rolling logical any, encoded as zero/one/`roll_na_logical`.
      integer, intent(in), optional :: min_obs !! Minimum nonmissing observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vectors have no cross-column row filtering.
      logical, intent(in), optional :: na_restore !! Restore a missing current observation when true.
      logical, intent(in), optional :: online !! Accepted for API parity; this implementation uses deterministic window evaluation.
      integer :: i, j, lo, n_obs, min_n
      logical :: restore, any_true, any_missing

      min_n = effective_min_obs(width, min_obs)
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x)))
      result = roll_na_logical
      do i = 1, size(x)
         if (restore .and. x(i) == roll_na_logical) cycle
         lo = max(1, i - width + 1)
         n_obs = 0
         any_true = .false.
         any_missing = .false.
         do j = lo, i
            if (x(j) == roll_na_logical) then
               any_missing = .true.
            else
               n_obs = n_obs + 1
               if (x(j) /= 0) any_true = .true.
            end if
         end do
         if (n_obs < min_n) cycle
         if (any_true) then
            result(i) = 1
         else if (any_missing) then
            result(i) = roll_na_logical
         else
            result(i) = 0
         end if
      end do
   end subroutine roll_any_vec

   pure subroutine roll_all_vec(x, width, result, min_obs, complete_obs, na_restore, online)
      integer, intent(in) :: x(:) !! Logical series encoded as zero/one with `roll_na_logical` for missing values.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      integer, allocatable, intent(out) :: result(:) !! Rolling logical all, encoded as zero/one/`roll_na_logical`.
      integer, intent(in), optional :: min_obs !! Minimum nonmissing observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vectors have no cross-column row filtering.
      logical, intent(in), optional :: na_restore !! Restore a missing current observation when true.
      logical, intent(in), optional :: online !! Accepted for API parity; this implementation uses deterministic window evaluation.
      integer :: i, j, lo, n_obs, min_n
      logical :: restore, any_false, any_missing

      min_n = effective_min_obs(width, min_obs)
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x)))
      result = roll_na_logical
      do i = 1, size(x)
         if (restore .and. x(i) == roll_na_logical) cycle
         lo = max(1, i - width + 1)
         n_obs = 0
         any_false = .false.
         any_missing = .false.
         do j = lo, i
            if (x(j) == roll_na_logical) then
               any_missing = .true.
            else
               n_obs = n_obs + 1
               if (x(j) == 0) any_false = .true.
            end if
         end do
         if (n_obs < min_n) cycle
         if (any_false) then
            result(i) = 0
         else if (any_missing) then
            result(i) = roll_na_logical
         else
            result(i) = 1
         end if
      end do
   end subroutine roll_all_vec

   pure subroutine roll_sum_vec(x, width, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:) !! Numeric time series; NaNs are treated as missing observations.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:) !! Weighted rolling sums.
      real(dp), intent(in), optional :: weights(:) !! Weights whose final element aligns with the current observation.
      integer, intent(in), optional :: min_obs !! Minimum finite observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vectors have no cross-column row filtering.
      logical, intent(in), optional :: na_restore !! Restore a NaN current observation when true.
      logical, intent(in), optional :: online !! Accepted for API parity; batch evaluation is used for both settings.
      integer :: i, c, n_obs, min_n
      real(dp) :: acc
      logical :: restore

      min_n = effective_min_obs(width, min_obs)
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x)))
      result = na_real()
      do i = 1, size(x)
         if (restore .and. ieee_is_nan(x(i))) cycle
         acc = 0.0_dp
         n_obs = 0
         do c = 0, min(width - 1, i - 1)
            if (.not. ieee_is_nan(x(i - c))) then
               acc = acc + window_weight(weights, c) * x(i - c)
               n_obs = n_obs + 1
            end if
         end do
         if (n_obs >= min_n) result(i) = acc
      end do
   end subroutine roll_sum_vec

   pure subroutine roll_prod_vec(x, width, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:) !! Numeric time series; NaNs are treated as missing observations.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:) !! Rolling products of weighted observations.
      real(dp), intent(in), optional :: weights(:) !! Weights whose final element aligns with the current observation.
      integer, intent(in), optional :: min_obs !! Minimum finite observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vectors have no cross-column row filtering.
      logical, intent(in), optional :: na_restore !! Restore a NaN current observation when true.
      logical, intent(in), optional :: online !! Accepted for API parity; batch evaluation is used for both settings.
      integer :: i, c, n_obs, min_n
      real(dp) :: acc
      logical :: restore

      min_n = effective_min_obs(width, min_obs)
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x)))
      result = na_real()
      do i = 1, size(x)
         if (restore .and. ieee_is_nan(x(i))) cycle
         acc = 1.0_dp
         n_obs = 0
         do c = 0, min(width - 1, i - 1)
            if (.not. ieee_is_nan(x(i - c))) then
               acc = acc * (window_weight(weights, c) * x(i - c))
               n_obs = n_obs + 1
            end if
         end do
         if (n_obs >= min_n) result(i) = acc
      end do
   end subroutine roll_prod_vec

   pure subroutine roll_mean_vec(x, width, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:) !! Numeric time series; NaNs are treated as missing observations.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:) !! Weighted rolling means.
      real(dp), intent(in), optional :: weights(:) !! Weights whose final element aligns with the current observation.
      integer, intent(in), optional :: min_obs !! Minimum finite observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vectors have no cross-column row filtering.
      logical, intent(in), optional :: na_restore !! Restore a NaN current observation when true.
      logical, intent(in), optional :: online !! Accepted for API parity; batch evaluation is used for both settings.
      integer :: i, c, n_obs, min_n
      real(dp) :: sum_w, sum_x, w
      logical :: restore

      min_n = effective_min_obs(width, min_obs)
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x)))
      result = na_real()
      do i = 1, size(x)
         if (restore .and. ieee_is_nan(x(i))) cycle
         sum_w = 0.0_dp
         sum_x = 0.0_dp
         n_obs = 0
         do c = 0, min(width - 1, i - 1)
            if (.not. ieee_is_nan(x(i - c))) then
               w = window_weight(weights, c)
               sum_w = sum_w + w
               sum_x = sum_x + w * x(i - c)
               n_obs = n_obs + 1
            end if
         end do
         if (n_obs >= min_n .and. sum_w /= 0.0_dp) result(i) = sum_x / sum_w
      end do
   end subroutine roll_mean_vec

   pure subroutine roll_min_vec(x, width, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:) !! Numeric time series; NaNs are treated as missing observations.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:) !! Rolling minima.
      real(dp), intent(in), optional :: weights(:) !! Accepted for R API parity; positive weights do not alter extrema.
      integer, intent(in), optional :: min_obs !! Minimum finite observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vectors have no cross-column row filtering.
      logical, intent(in), optional :: na_restore !! Restore a NaN current observation when true.
      logical, intent(in), optional :: online !! Accepted for API parity; batch evaluation is used for both settings.
      call roll_extreme_vec(x, width, result, .true., min_obs, na_restore)
   end subroutine roll_min_vec

   pure subroutine roll_max_vec(x, width, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:) !! Numeric time series; NaNs are treated as missing observations.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:) !! Rolling maxima.
      real(dp), intent(in), optional :: weights(:) !! Accepted for R API parity; positive weights do not alter extrema.
      integer, intent(in), optional :: min_obs !! Minimum finite observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vectors have no cross-column row filtering.
      logical, intent(in), optional :: na_restore !! Restore a NaN current observation when true.
      logical, intent(in), optional :: online !! Accepted for API parity; batch evaluation is used for both settings.
      call roll_extreme_vec(x, width, result, .false., min_obs, na_restore)
   end subroutine roll_max_vec

   pure subroutine roll_extreme_vec(x, width, result, want_min, min_obs, na_restore)
      real(dp), intent(in) :: x(:) !! Numeric series to scan.
      integer, intent(in) :: width !! Rolling window width.
      real(dp), allocatable, intent(out) :: result(:) !! Rolling extrema.
      logical, intent(in) :: want_min !! Select minimum when true and maximum when false.
      integer, intent(in), optional :: min_obs !! Minimum finite observations.
      logical, intent(in), optional :: na_restore !! Restore a NaN current observation when true.
      integer :: i, c, n_obs, min_n
      real(dp) :: value
      logical :: restore, first

      min_n = effective_min_obs(width, min_obs)
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x)))
      result = na_real()
      do i = 1, size(x)
         if (restore .and. ieee_is_nan(x(i))) cycle
         n_obs = 0
         first = .true.
         value = 0.0_dp
         do c = 0, min(width - 1, i - 1)
            if (.not. ieee_is_nan(x(i - c))) then
               if (first) then
                  value = x(i - c)
                  first = .false.
               else if (want_min) then
                  value = min(value, x(i - c))
               else
                  value = max(value, x(i - c))
               end if
               n_obs = n_obs + 1
            end if
         end do
         if (n_obs >= min_n) result(i) = value
      end do
   end subroutine roll_extreme_vec

   pure subroutine roll_idxmin_vec(x, width, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:) !! Numeric time series; NaNs are treated as missing observations.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      integer, allocatable, intent(out) :: result(:) !! One-based index within each current window of its minimum.
      real(dp), intent(in), optional :: weights(:) !! Accepted for R API parity; positive weights do not alter the index.
      integer, intent(in), optional :: min_obs !! Minimum finite observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vectors have no cross-column row filtering.
      logical, intent(in), optional :: na_restore !! Return `roll_na_logical` when the current value is missing and true.
      logical, intent(in), optional :: online !! Accepted for API parity; batch evaluation is used for both settings.
      call roll_index_extreme_vec(x, width, result, .true., min_obs, na_restore)
   end subroutine roll_idxmin_vec

   pure subroutine roll_idxmax_vec(x, width, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:) !! Numeric time series; NaNs are treated as missing observations.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      integer, allocatable, intent(out) :: result(:) !! One-based index within each current window of its maximum.
      real(dp), intent(in), optional :: weights(:) !! Accepted for R API parity; positive weights do not alter the index.
      integer, intent(in), optional :: min_obs !! Minimum finite observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vectors have no cross-column row filtering.
      logical, intent(in), optional :: na_restore !! Return `roll_na_logical` when the current value is missing and true.
      logical, intent(in), optional :: online !! Accepted for API parity; batch evaluation is used for both settings.
      call roll_index_extreme_vec(x, width, result, .false., min_obs, na_restore)
   end subroutine roll_idxmax_vec

   pure subroutine roll_index_extreme_vec(x, width, result, want_min, min_obs, na_restore)
      real(dp), intent(in) :: x(:) !! Numeric series to scan.
      integer, intent(in) :: width !! Rolling window width.
      integer, allocatable, intent(out) :: result(:) !! Window-relative one-based extrema indices.
      logical, intent(in) :: want_min !! Select minimum when true and maximum when false.
      integer, intent(in), optional :: min_obs !! Minimum finite observations.
      logical, intent(in), optional :: na_restore !! Restore missing current rows when true.
      integer :: i, c, idx, lo, n_obs, min_n
      real(dp) :: value
      logical :: restore, first

      min_n = effective_min_obs(width, min_obs)
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x)))
      result = roll_na_logical
      do i = 1, size(x)
         if (restore .and. ieee_is_nan(x(i))) cycle
         lo = max(1, i - width + 1)
         first = .true.
         n_obs = 0
         idx = i
         value = 0.0_dp
         do c = 0, min(width - 1, i - 1)
            if (.not. ieee_is_nan(x(i - c))) then
               if (first) then
                  value = x(i - c)
                  idx = i - c
                  first = .false.
               else if (want_min) then
                  if (x(i - c) <= value) then
                     value = x(i - c)
                     idx = i - c
                  end if
               else
                  if (x(i - c) >= value) then
                     value = x(i - c)
                     idx = i - c
                  end if
               end if
               n_obs = n_obs + 1
            end if
         end do
         if (n_obs >= min_n) result(i) = idx - lo + 1
      end do
   end subroutine roll_index_extreme_vec

   pure subroutine roll_median_vec(x, width, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:) !! Numeric time series; NaNs are treated as missing observations.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:) !! Rolling weighted medians using Hyndman-Fan type 2 semantics.
      real(dp), intent(in), optional :: weights(:) !! Positive observation weights aligned from the end of the vector.
      integer, intent(in), optional :: min_obs !! Minimum finite observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vectors have no cross-column row filtering.
      logical, intent(in), optional :: na_restore !! Restore a NaN current observation when true.
      logical, intent(in), optional :: online !! Accepted for API parity; weighted batch evaluation is used.
      call roll_quantile_vec(x, width, 0.5_dp, result, weights, min_obs, complete_obs, na_restore, online)
   end subroutine roll_median_vec

   pure subroutine roll_quantile_vec(x, width, p, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:) !! Numeric time series; NaNs are treated as missing observations.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), intent(in) :: p !! Quantile probability in the closed interval `[0,1]`.
      real(dp), allocatable, intent(out) :: result(:) !! Rolling weighted type-2 quantiles.
      real(dp), intent(in), optional :: weights(:) !! Positive observation weights aligned from the end of the vector.
      integer, intent(in), optional :: min_obs !! Minimum finite observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vectors have no cross-column row filtering.
      logical, intent(in), optional :: na_restore !! Restore a NaN current observation when true.
      logical, intent(in), optional :: online !! Accepted for API parity; weighted batch evaluation is used.
      integer :: i, c, j, k, n_obs, min_n
      real(dp) :: total_w, cum_w, temp_v, temp_w, target, tol
      real(dp), allocatable :: values(:), wts(:)
      logical :: restore

      min_n = effective_min_obs(width, min_obs)
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x)), values(width), wts(width))
      result = na_real()
      tol = sqrt(epsilon(1.0_dp))
      do i = 1, size(x)
         if (restore .and. ieee_is_nan(x(i))) cycle
         n_obs = 0
         do c = min(width - 1, i - 1), 0, -1
            if (.not. ieee_is_nan(x(i - c))) then
               n_obs = n_obs + 1
               values(n_obs) = x(i - c)
               wts(n_obs) = window_weight(weights, c)
            end if
         end do
         if (n_obs < min_n) cycle
         do j = 2, n_obs
            temp_v = values(j)
            temp_w = wts(j)
            k = j - 1
            do while (k >= 1)
               if (values(k) <= temp_v) exit
               values(k + 1) = values(k)
               wts(k + 1) = wts(k)
               k = k - 1
            end do
            values(k + 1) = temp_v
            wts(k + 1) = temp_w
         end do
         total_w = sum(wts(1:n_obs))
         if (total_w <= 0.0_dp) cycle
         if (p <= 0.0_dp) then
            result(i) = values(1)
            cycle
         else if (p >= 1.0_dp) then
            result(i) = values(n_obs)
            cycle
         end if
         target = p * total_w
         cum_w = 0.0_dp
         do j = 1, n_obs
            cum_w = cum_w + wts(j)
            if (cum_w >= target) then
               if (abs(cum_w / total_w - p) <= tol .and. j < n_obs) then
                  result(i) = 0.5_dp * (values(j) + values(j + 1))
               else
                  result(i) = values(j)
               end if
               exit
            end if
         end do
      end do
   end subroutine roll_quantile_vec

   pure subroutine roll_var_vec(x, width, result, weights, center, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:) !! Numeric time series; NaNs are treated as missing observations.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:) !! Unbiased weighted rolling variances.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end of the vector.
      logical, intent(in), optional :: center !! Use the weighted mean when true; use zero when false.
      integer, intent(in), optional :: min_obs !! Minimum finite observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vectors have no cross-column row filtering.
      logical, intent(in), optional :: na_restore !! Restore a NaN current observation when true.
      logical, intent(in), optional :: online !! Accepted for API parity; batch evaluation is used.
      integer :: i, c, n_obs, min_n
      real(dp) :: sum_w, sum_w2, sum_x, mean_x, ss, w, denom, value
      logical :: use_center, restore

      min_n = effective_min_obs(width, min_obs)
      use_center = .true.
      if (present(center)) use_center = center
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x)))
      result = na_real()
      do i = 1, size(x)
         if (restore .and. ieee_is_nan(x(i))) cycle
         sum_w = 0.0_dp
         sum_w2 = 0.0_dp
         sum_x = 0.0_dp
         n_obs = 0
         do c = 0, min(width - 1, i - 1)
            if (.not. ieee_is_nan(x(i - c))) then
               w = window_weight(weights, c)
               sum_w = sum_w + w
               sum_w2 = sum_w2 + w * w
               sum_x = sum_x + w * x(i - c)
               n_obs = n_obs + 1
            end if
         end do
         if (n_obs <= 1 .or. n_obs < min_n .or. sum_w == 0.0_dp) cycle
         mean_x = 0.0_dp
         if (use_center) mean_x = sum_x / sum_w
         ss = 0.0_dp
         do c = 0, min(width - 1, i - 1)
            if (.not. ieee_is_nan(x(i - c))) then
               w = window_weight(weights, c)
               ss = ss + w * (x(i - c) - mean_x) ** 2
            end if
         end do
         denom = sum_w - sum_w2 / sum_w
         if (denom == 0.0_dp) cycle
         value = ss / denom
         if (value > epsilon(1.0_dp)) then
            result(i) = value
         else if (value > -epsilon(1.0_dp)) then
            result(i) = 0.0_dp
         end if
      end do
   end subroutine roll_var_vec

   pure subroutine roll_sd_vec(x, width, result, weights, center, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:) !! Numeric time series; NaNs are treated as missing observations.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:) !! Rolling standard deviations.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end of the vector.
      logical, intent(in), optional :: center !! Use the weighted mean when true; use zero when false.
      integer, intent(in), optional :: min_obs !! Minimum finite observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vectors have no cross-column row filtering.
      logical, intent(in), optional :: na_restore !! Restore a NaN current observation when true.
      logical, intent(in), optional :: online !! Accepted for API parity; batch evaluation is used.
      real(dp), allocatable :: variance(:)
      integer :: i

      call roll_var_vec(x, width, variance, weights, center, min_obs, complete_obs, na_restore, online)
      allocate(result(size(x)))
      result = na_real()
      do i = 1, size(x)
         if (.not. ieee_is_nan(variance(i))) result(i) = sqrt(variance(i))
      end do
   end subroutine roll_sd_vec

   pure subroutine roll_scale_vec(x, width, result, weights, center, scale, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:) !! Numeric time series; NaNs are treated as missing observations.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:) !! Latest finite observation after rolling centering/scaling.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end of the vector.
      logical, intent(in), optional :: center !! Subtract the weighted rolling mean when true.
      logical, intent(in), optional :: scale !! Divide by weighted rolling standard deviation when true.
      integer, intent(in), optional :: min_obs !! Minimum finite observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Accepted for R API parity; vectors have no cross-column row filtering.
      logical, intent(in), optional :: na_restore !! Restore a NaN current observation when true.
      logical, intent(in), optional :: online !! Accepted for API parity; batch evaluation is used.
      integer :: i, c, n_obs, min_n
      real(dp) :: sum_w, sum_w2, sum_x, mean_x, ss, variance, w, latest
      logical :: use_center, use_scale, restore, found

      min_n = effective_min_obs(width, min_obs)
      use_center = .true.
      if (present(center)) use_center = center
      use_scale = .true.
      if (present(scale)) use_scale = scale
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x)))
      result = na_real()
      do i = 1, size(x)
         if (restore .and. ieee_is_nan(x(i))) cycle
         sum_w = 0.0_dp
         sum_w2 = 0.0_dp
         sum_x = 0.0_dp
         n_obs = 0
         found = .false.
         latest = 0.0_dp
         do c = 0, min(width - 1, i - 1)
            if (.not. ieee_is_nan(x(i - c))) then
               if (.not. found) then
                  latest = x(i - c)
                  found = .true.
               end if
               w = window_weight(weights, c)
               sum_w = sum_w + w
               sum_w2 = sum_w2 + w * w
               sum_x = sum_x + w * x(i - c)
               n_obs = n_obs + 1
            end if
         end do
         if (n_obs < min_n .or. .not. found .or. sum_w == 0.0_dp) cycle
         mean_x = 0.0_dp
         if (use_center) mean_x = sum_x / sum_w
         if (.not. use_scale) then
            result(i) = latest - mean_x
            cycle
         end if
         if (n_obs <= 1) cycle
         ss = 0.0_dp
         do c = 0, min(width - 1, i - 1)
            if (.not. ieee_is_nan(x(i - c))) then
               w = window_weight(weights, c)
               ss = ss + w * (x(i - c) - mean_x) ** 2
            end if
         end do
         variance = ss / (sum_w - sum_w2 / sum_w)
         if (variance > epsilon(1.0_dp)) result(i) = (latest - mean_x) / sqrt(variance)
      end do
   end subroutine roll_scale_vec

   pure subroutine roll_pair_vec(x, y, width, result, weights, center, scale, min_obs, na_restore, crossprod_mode)
      real(dp), intent(in) :: x(:) !! First numeric series; NaNs are pairwise-missing.
      real(dp), intent(in) :: y(:) !! Second numeric series with the same length as `x`.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:) !! Rolling covariance/correlation/crossproduct values.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end of the vector.
      logical, intent(in), optional :: center !! Subtract weighted means when true.
      logical, intent(in), optional :: scale !! Normalize to correlation when true.
      integer, intent(in), optional :: min_obs !! Minimum pairwise-complete observations; defaults to `width`.
      logical, intent(in), optional :: na_restore !! Restore a NaN at the current pair when true.
      logical, intent(in), optional :: crossprod_mode !! Omit unbiased covariance divisor when true.
      integer :: i, c, n_obs, min_n
      real(dp) :: sum_w, sum_w2, sx, sy, mx, my, ssx, ssy, sxy, w, denom
      logical :: use_center, use_scale, restore, cross_mode

      min_n = effective_min_obs(width, min_obs)
      use_center = .true.
      if (present(center)) use_center = center
      use_scale = .false.
      if (present(scale)) use_scale = scale
      restore = .false.
      if (present(na_restore)) restore = na_restore
      cross_mode = .false.
      if (present(crossprod_mode)) cross_mode = crossprod_mode
      allocate(result(size(x)))
      result = na_real()
      do i = 1, size(x)
         if (restore .and. (ieee_is_nan(x(i)) .or. ieee_is_nan(y(i)))) cycle
         sum_w = 0.0_dp
         sum_w2 = 0.0_dp
         sx = 0.0_dp
         sy = 0.0_dp
         n_obs = 0
         do c = 0, min(width - 1, i - 1)
            if (.not. ieee_is_nan(x(i - c)) .and. .not. ieee_is_nan(y(i - c))) then
               w = window_weight(weights, c)
               sum_w = sum_w + w
               sum_w2 = sum_w2 + w * w
               sx = sx + w * x(i - c)
               sy = sy + w * y(i - c)
               n_obs = n_obs + 1
            end if
         end do
         if (n_obs < min_n .or. sum_w == 0.0_dp) cycle
         if (.not. cross_mode .and. n_obs <= 1) cycle
         mx = 0.0_dp
         my = 0.0_dp
         if (use_center) then
            mx = sx / sum_w
            my = sy / sum_w
         end if
         ssx = 0.0_dp
         ssy = 0.0_dp
         sxy = 0.0_dp
         do c = 0, min(width - 1, i - 1)
            if (.not. ieee_is_nan(x(i - c)) .and. .not. ieee_is_nan(y(i - c))) then
               w = window_weight(weights, c)
               ssx = ssx + w * (x(i - c) - mx) ** 2
               ssy = ssy + w * (y(i - c) - my) ** 2
               sxy = sxy + w * (x(i - c) - mx) * (y(i - c) - my)
            end if
         end do
         if (use_scale) then
            if (ssx > epsilon(1.0_dp) .and. ssy > epsilon(1.0_dp)) result(i) = sxy / sqrt(ssx * ssy)
         else if (cross_mode) then
            result(i) = sxy
         else
            denom = sum_w - sum_w2 / sum_w
            if (denom /= 0.0_dp) result(i) = sxy / denom
         end if
      end do
   end subroutine roll_pair_vec

end module roll_univariate
