! SPDX-License-Identifier: MIT
module r_descriptive
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_negative_inf, &
      ieee_positive_inf, ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   implicit none
   private

   public :: r_correlation, r_count_nonmissing, r_covariance
   public :: r_mean, r_variance, r_sd
   public :: r_weighted_correlation, r_weighted_covariance, r_weighted_mean
   public :: r_weighted_sd, r_weighted_variance

contains

   pure integer function r_count_nonmissing(x) result(n)
      real(dp), intent(in) :: x(:)
      integer :: i

      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) n = n + 1
      end do
   end function r_count_nonmissing

   pure real(dp) function r_mean(x, na_rm) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: na_rm
      real(dp) :: total, correction, y, updated
      integer :: i, n
      logical :: remove_na, has_positive_inf, has_negative_inf

      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      total = 0.0_dp
      correction = 0.0_dp
      n = 0
      has_positive_inf = .false.
      has_negative_inf = .false.

      do i = 1, size(x)
         if (ieee_is_nan(x(i))) then
            if (remove_na) cycle
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
         n = n + 1
         if (.not. ieee_is_finite(x(i))) then
            if (x(i) > 0.0_dp) has_positive_inf = .true.
            if (x(i) < 0.0_dp) has_negative_inf = .true.
            cycle
         end if
         y = x(i) - correction
         updated = total + y
         correction = (updated - total) - y
         total = updated
      end do

      if (n == 0 .or. (has_positive_inf .and. has_negative_inf)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (has_positive_inf) then
         value = ieee_value(0.0_dp, ieee_positive_inf)
      else if (has_negative_inf) then
         value = ieee_value(0.0_dp, ieee_negative_inf)
      else
         value = total / real(n, dp)
      end if
   end function r_mean

   pure real(dp) function r_variance(x, na_rm, center, ddof) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: na_rm
      real(dp), intent(in), optional :: center
      integer, intent(in), optional :: ddof
      real(dp) :: mean_value, total, correction, term, updated
      integer :: degrees, i, n
      logical :: remove_na

      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      degrees = 1
      if (present(ddof)) degrees = ddof
      if (degrees < 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      if (present(center)) then
         mean_value = center
      else
         mean_value = r_mean(x, remove_na)
      end if
      if (.not. ieee_is_finite(mean_value)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      total = 0.0_dp
      correction = 0.0_dp
      n = 0
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) then
            if (remove_na) cycle
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
         if (.not. ieee_is_finite(x(i))) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
         n = n + 1
         term = (x(i) - mean_value)**2 - correction
         updated = total + term
         correction = (updated - total) - term
         total = updated
      end do

      if (n <= degrees) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         value = total / real(n - degrees, dp)
      end if
   end function r_variance

   pure real(dp) function r_sd(x, na_rm, center, ddof) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: na_rm
      real(dp), intent(in), optional :: center
      integer, intent(in), optional :: ddof

      value = r_variance(x, na_rm, center, ddof)
      if (value >= 0.0_dp) value = sqrt(value)
   end function r_sd

   pure real(dp) function r_covariance(x, y, na_rm, finite_only, ddof) result(value)
      real(dp), intent(in) :: x(:), y(:)
      logical, intent(in), optional :: na_rm, finite_only
      integer, intent(in), optional :: ddof
      real(dp) :: mean_x, mean_y, total_x, total_y
      integer :: degrees, i, n
      logical :: remove_na, require_finite

      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(x) /= size(y)) return
      degrees = 1
      if (present(ddof)) degrees = ddof
      if (degrees < 0) return
      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      require_finite = .false.
      if (present(finite_only)) require_finite = finite_only

      n = 0
      total_x = 0.0_dp
      total_y = 0.0_dp
      do i = 1, size(x)
         if (require_finite .and. (.not. ieee_is_finite(x(i)) .or. .not. ieee_is_finite(y(i)))) cycle
         if (remove_na .and. (ieee_is_nan(x(i)) .or. ieee_is_nan(y(i)))) cycle
         if (.not. ieee_is_finite(x(i)) .or. .not. ieee_is_finite(y(i))) return
         n = n + 1
         total_x = total_x + x(i)
         total_y = total_y + y(i)
      end do
      if (n <= degrees) return
      mean_x = total_x/real(n, dp)
      mean_y = total_y/real(n, dp)

      value = 0.0_dp
      do i = 1, size(x)
         if (require_finite .and. (.not. ieee_is_finite(x(i)) .or. .not. ieee_is_finite(y(i)))) cycle
         if (remove_na .and. (ieee_is_nan(x(i)) .or. ieee_is_nan(y(i)))) cycle
         value = value + (x(i) - mean_x)*(y(i) - mean_y)
      end do
      value = value/real(n - degrees, dp)
   end function r_covariance

   pure real(dp) function r_correlation(x, y, na_rm, finite_only) result(value)
      real(dp), intent(in) :: x(:), y(:)
      logical, intent(in), optional :: na_rm, finite_only
      real(dp) :: mean_x, mean_y, sum_xx, sum_xy, sum_yy, total_x, total_y
      integer :: i, n
      logical :: remove_na, require_finite

      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(x) /= size(y)) return
      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      require_finite = .false.
      if (present(finite_only)) require_finite = finite_only

      n = 0
      total_x = 0.0_dp
      total_y = 0.0_dp
      do i = 1, size(x)
         if (require_finite .and. (.not. ieee_is_finite(x(i)) .or. .not. ieee_is_finite(y(i)))) cycle
         if (remove_na .and. (ieee_is_nan(x(i)) .or. ieee_is_nan(y(i)))) cycle
         if (.not. ieee_is_finite(x(i)) .or. .not. ieee_is_finite(y(i))) return
         n = n + 1
         total_x = total_x + x(i)
         total_y = total_y + y(i)
      end do
      if (n < 2) return
      mean_x = total_x/real(n, dp)
      mean_y = total_y/real(n, dp)

      sum_xx = 0.0_dp
      sum_xy = 0.0_dp
      sum_yy = 0.0_dp
      do i = 1, size(x)
         if (require_finite .and. (.not. ieee_is_finite(x(i)) .or. .not. ieee_is_finite(y(i)))) cycle
         if (remove_na .and. (ieee_is_nan(x(i)) .or. ieee_is_nan(y(i)))) cycle
         sum_xx = sum_xx + (x(i) - mean_x)**2
         sum_xy = sum_xy + (x(i) - mean_x)*(y(i) - mean_y)
         sum_yy = sum_yy + (y(i) - mean_y)**2
      end do
      if (sum_xx <= 0.0_dp .or. sum_yy <= 0.0_dp) return
      value = sum_xy/sqrt(sum_xx*sum_yy)
   end function r_correlation

   pure real(dp) function r_weighted_mean(x, weights, na_rm, finite_only) result(value)
      real(dp), intent(in) :: x(:), weights(:)
      logical, intent(in), optional :: na_rm, finite_only
      real(dp) :: total, total_weight
      integer :: i
      logical :: remove_na, require_finite

      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. valid_weights(weights, size(x))) return
      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      require_finite = .false.
      if (present(finite_only)) require_finite = finite_only
      total = 0.0_dp
      total_weight = 0.0_dp
      do i = 1, size(x)
         if (weights(i) == 0.0_dp) cycle
         if (require_finite .and. .not. ieee_is_finite(x(i))) cycle
         if (remove_na .and. ieee_is_nan(x(i))) cycle
         if (ieee_is_nan(x(i))) return
         total = total + weights(i)*x(i)
         total_weight = total_weight + weights(i)
      end do
      if (total_weight <= 0.0_dp) return
      value = total/total_weight
   end function r_weighted_mean

   pure real(dp) function r_weighted_variance(x, weights, unbiased, na_rm, finite_only) result(value)
      real(dp), intent(in) :: x(:), weights(:)
      logical, intent(in), optional :: unbiased, na_rm, finite_only
      real(dp) :: denominator, mean_value, sum_weight_squared, total_weight
      integer :: i
      logical :: correct_bias, remove_na, require_finite

      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. valid_weights(weights, size(x))) return
      correct_bias = .false.
      if (present(unbiased)) correct_bias = unbiased
      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      require_finite = .false.
      if (present(finite_only)) require_finite = finite_only
      mean_value = r_weighted_mean(x, weights, remove_na, require_finite)
      if (.not. ieee_is_finite(mean_value)) return

      total_weight = 0.0_dp
      sum_weight_squared = 0.0_dp
      value = 0.0_dp
      do i = 1, size(x)
         if (weights(i) == 0.0_dp) cycle
         if (require_finite .and. .not. ieee_is_finite(x(i))) cycle
         if (remove_na .and. ieee_is_nan(x(i))) cycle
         total_weight = total_weight + weights(i)
         sum_weight_squared = sum_weight_squared + weights(i)**2
         value = value + weights(i)*(x(i) - mean_value)**2
      end do
      denominator = total_weight
      if (correct_bias) denominator = total_weight - sum_weight_squared/total_weight
      if (denominator <= 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         value = value/denominator
      end if
   end function r_weighted_variance

   pure real(dp) function r_weighted_sd(x, weights, unbiased, na_rm, finite_only) result(value)
      real(dp), intent(in) :: x(:), weights(:)
      logical, intent(in), optional :: unbiased, na_rm, finite_only

      value = r_weighted_variance(x, weights, unbiased, na_rm, finite_only)
      if (value >= 0.0_dp) value = sqrt(value)
   end function r_weighted_sd

   pure real(dp) function r_weighted_covariance(x, y, weights, unbiased, na_rm, finite_only) result(value)
      real(dp), intent(in) :: x(:), y(:), weights(:)
      logical, intent(in), optional :: unbiased, na_rm, finite_only
      real(dp) :: denominator, mean_x, mean_y, sum_weight_squared, total_weight
      integer :: i
      logical :: correct_bias, remove_na, require_finite

      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(x) /= size(y) .or. .not. valid_weights(weights, size(x))) return
      correct_bias = .false.
      if (present(unbiased)) correct_bias = unbiased
      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      require_finite = .false.
      if (present(finite_only)) require_finite = finite_only
      call weighted_pair_means(x, y, weights, remove_na, require_finite, mean_x, mean_y, total_weight)
      if (.not. ieee_is_finite(mean_x) .or. .not. ieee_is_finite(mean_y)) return

      sum_weight_squared = 0.0_dp
      value = 0.0_dp
      do i = 1, size(x)
         if (.not. keep_weighted_pair(x(i), y(i), weights(i), remove_na, require_finite)) cycle
         sum_weight_squared = sum_weight_squared + weights(i)**2
         value = value + weights(i)*(x(i) - mean_x)*(y(i) - mean_y)
      end do
      denominator = total_weight
      if (correct_bias) denominator = total_weight - sum_weight_squared/total_weight
      if (denominator <= 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         value = value/denominator
      end if
   end function r_weighted_covariance

   pure real(dp) function r_weighted_correlation(x, y, weights, na_rm, finite_only) result(value)
      real(dp), intent(in) :: x(:), y(:), weights(:)
      logical, intent(in), optional :: na_rm, finite_only
      real(dp) :: mean_x, mean_y, sum_xx, sum_xy, sum_yy, total_weight
      integer :: i
      logical :: remove_na, require_finite

      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(x) /= size(y) .or. .not. valid_weights(weights, size(x))) return
      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      require_finite = .false.
      if (present(finite_only)) require_finite = finite_only
      call weighted_pair_means(x, y, weights, remove_na, require_finite, mean_x, mean_y, total_weight)
      if (.not. ieee_is_finite(mean_x) .or. .not. ieee_is_finite(mean_y)) return

      sum_xx = 0.0_dp
      sum_xy = 0.0_dp
      sum_yy = 0.0_dp
      do i = 1, size(x)
         if (.not. keep_weighted_pair(x(i), y(i), weights(i), remove_na, require_finite)) cycle
         sum_xx = sum_xx + weights(i)*(x(i) - mean_x)**2
         sum_xy = sum_xy + weights(i)*(x(i) - mean_x)*(y(i) - mean_y)
         sum_yy = sum_yy + weights(i)*(y(i) - mean_y)**2
      end do
      if (sum_xx <= 0.0_dp .or. sum_yy <= 0.0_dp) return
      value = sum_xy/sqrt(sum_xx*sum_yy)
   end function r_weighted_correlation

   pure logical function valid_weights(weights, expected_size) result(valid)
      real(dp), intent(in) :: weights(:)
      integer, intent(in) :: expected_size

      valid = size(weights) == expected_size
      if (.not. valid) return
      valid = all(ieee_is_finite(weights)) .and. all(weights >= 0.0_dp)
   end function valid_weights

   pure logical function keep_weighted_pair(x, y, weight, remove_na, require_finite) result(keep)
      real(dp), intent(in) :: x, y, weight
      logical, intent(in) :: remove_na, require_finite

      keep = weight > 0.0_dp
      if (.not. keep) return
      if (require_finite .and. (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(y))) keep = .false.
      if (remove_na .and. (ieee_is_nan(x) .or. ieee_is_nan(y))) keep = .false.
   end function keep_weighted_pair

   pure subroutine weighted_pair_means(x, y, weights, remove_na, require_finite, mean_x, mean_y, total_weight)
      real(dp), intent(in) :: x(:), y(:), weights(:)
      logical, intent(in) :: remove_na, require_finite
      real(dp), intent(out) :: mean_x, mean_y, total_weight
      integer :: i

      mean_x = 0.0_dp
      mean_y = 0.0_dp
      total_weight = 0.0_dp
      do i = 1, size(x)
         if (.not. keep_weighted_pair(x(i), y(i), weights(i), remove_na, require_finite)) cycle
         if (ieee_is_nan(x(i)) .or. ieee_is_nan(y(i))) then
            mean_x = ieee_value(0.0_dp, ieee_quiet_nan)
            mean_y = mean_x
            return
         end if
         mean_x = mean_x + weights(i)*x(i)
         mean_y = mean_y + weights(i)*y(i)
         total_weight = total_weight + weights(i)
      end do
      if (total_weight <= 0.0_dp) then
         mean_x = ieee_value(0.0_dp, ieee_quiet_nan)
         mean_y = mean_x
      else
         mean_x = mean_x/total_weight
         mean_y = mean_y/total_weight
      end if
   end subroutine weighted_pair_means

end module r_descriptive
