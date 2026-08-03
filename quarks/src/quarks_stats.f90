module quarks_stats
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use quarks_kinds, only : dp
   implicit none
   private

   public :: mean_value, sample_variance, quantile_type7, weighted_quantile
   public :: sort_real, sort_real_with_index, chi_square_upper_tail
   public :: binomial_cdf, nan_value, safe_log_probability

contains

   pure function nan_value() result(value)
      real(dp) :: value
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_value

   pure function mean_value(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      if (size(x) == 0) then
         value = nan_value()
      else
         value = sum(x) / real(size(x), dp)
      end if
   end function mean_value

   pure function sample_variance(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value, center
      if (size(x) <= 1) then
         value = nan_value()
         return
      end if
      center = mean_value(x)
      value = sum((x - center)**2) / real(size(x) - 1, dp)
   end function sample_variance

   recursive subroutine quicksort_pairs(values, indices, left, right)
      real(dp), intent(inout) :: values(:)
      integer, intent(inout) :: indices(:)
      integer, intent(in) :: left, right
      integer :: i, j, itmp
      real(dp) :: pivot, tmp
      if (left >= right) return
      i = left
      j = right
      pivot = values((left + right) / 2)
      do
         do while (values(i) < pivot)
            i = i + 1
         end do
         do while (values(j) > pivot)
            j = j - 1
         end do
         if (i <= j) then
            tmp = values(i)
            values(i) = values(j)
            values(j) = tmp
            itmp = indices(i)
            indices(i) = indices(j)
            indices(j) = itmp
            i = i + 1
            j = j - 1
         end if
         if (i > j) exit
      end do
      if (left < j) call quicksort_pairs(values, indices, left, j)
      if (i < right) call quicksort_pairs(values, indices, i, right)
   end subroutine quicksort_pairs

   subroutine sort_real_with_index(x, sorted, indices)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: sorted(:)
      integer, allocatable, intent(out) :: indices(:)
      integer :: i
      allocate(sorted(size(x)), indices(size(x)))
      sorted = x
      do i = 1, size(x)
         indices(i) = i
      end do
      if (size(x) > 1) call quicksort_pairs(sorted, indices, 1, size(x))
   end subroutine sort_real_with_index

   subroutine sort_real(x, sorted)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: sorted(:)
      integer, allocatable :: indices(:)
      call sort_real_with_index(x, sorted, indices)
   end subroutine sort_real

   function quantile_type7(x, p) result(value)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: p
      real(dp) :: value, h, gamma
      real(dp), allocatable :: sorted(:)
      integer :: lo, hi, n
      n = size(x)
      if (n == 0 .or. p < 0.0_dp .or. p > 1.0_dp) then
         value = nan_value()
         return
      end if
      call sort_real(x, sorted)
      if (n == 1 .or. p <= 0.0_dp) then
         value = sorted(1)
         return
      end if
      if (p >= 1.0_dp) then
         value = sorted(n)
         return
      end if
      h = 1.0_dp + real(n - 1, dp) * p
      lo = int(floor(h))
      hi = int(ceiling(h))
      gamma = h - real(lo, dp)
      value = (1.0_dp - gamma) * sorted(lo) + gamma * sorted(hi)
   end function quantile_type7

   function weighted_quantile(x, weights, p) result(value)
      real(dp), intent(in) :: x(:), weights(:)
      real(dp), intent(in) :: p
      real(dp) :: value, total, target, fraction
      real(dp), allocatable :: sorted(:), sorted_weights(:), cumulative(:)
      integer, allocatable :: indices(:)
      integer :: i, high
      if (size(x) == 0 .or. size(x) /= size(weights) .or. &
          p < 0.0_dp .or. p > 1.0_dp .or. any(weights < 0.0_dp)) then
         value = nan_value()
         return
      end if
      call sort_real_with_index(x, sorted, indices)
      allocate(sorted_weights(size(x)), cumulative(size(x)))
      do i = 1, size(x)
         sorted_weights(i) = weights(indices(i))
      end do
      total = sum(sorted_weights)
      if (total <= 0.0_dp) then
         value = nan_value()
         return
      end if
      sorted_weights = sorted_weights / total
      cumulative(1) = sorted_weights(1)
      do i = 2, size(x)
         cumulative(i) = cumulative(i - 1) + sorted_weights(i)
      end do
      target = p
      if (target <= cumulative(1)) then
         value = sorted(1)
         return
      end if
      high = size(x)
      do i = 2, size(x)
         if (cumulative(i) > target) then
            high = i
            exit
         end if
      end do
      if (high == size(x) .and. target >= cumulative(high)) then
         value = sorted(high)
         return
      end if
      fraction = (target - cumulative(high - 1)) / &
         max(cumulative(high) - cumulative(high - 1), tiny(1.0_dp))
      value = sorted(high - 1) + fraction * (sorted(high) - sorted(high - 1))
   end function weighted_quantile

   pure function chi_square_upper_tail(x, degrees_freedom) result(value)
      real(dp), intent(in) :: x
      integer, intent(in) :: degrees_freedom
      real(dp) :: value
      if (x <= 0.0_dp) then
         value = 1.0_dp
      else if (degrees_freedom == 1) then
         value = erfc(sqrt(0.5_dp * x))
      else if (degrees_freedom == 2) then
         value = exp(-0.5_dp * x)
      else
         value = nan_value()
      end if
   end function chi_square_upper_tail

   pure function safe_log_probability(probability, count) result(value)
      real(dp), intent(in) :: probability
      integer, intent(in) :: count
      real(dp) :: value
      if (count == 0) then
         value = 0.0_dp
      else if (probability <= 0.0_dp) then
         value = -huge(1.0_dp)
      else
         value = real(count, dp) * log(probability)
      end if
   end function safe_log_probability

   function binomial_cdf(k, n, probability) result(value)
      integer, intent(in) :: k, n
      real(dp), intent(in) :: probability
      real(dp) :: value, maximum_log, term_log
      real(dp), allocatable :: logs(:)
      integer :: j, upper
      if (n < 0 .or. probability < 0.0_dp .or. probability > 1.0_dp) then
         value = nan_value()
         return
      end if
      if (k < 0) then
         value = 0.0_dp
         return
      end if
      if (k >= n) then
         value = 1.0_dp
         return
      end if
      if (probability <= 0.0_dp) then
         value = 1.0_dp
         return
      end if
      if (probability >= 1.0_dp) then
         value = 0.0_dp
         return
      end if
      upper = min(k, n)
      allocate(logs(0:upper))
      maximum_log = -huge(1.0_dp)
      do j = 0, upper
         term_log = log_gamma(real(n + 1, dp)) - log_gamma(real(j + 1, dp)) - &
            log_gamma(real(n - j + 1, dp)) + real(j, dp) * log(probability) + &
            real(n - j, dp) * log(1.0_dp - probability)
         logs(j) = term_log
         maximum_log = max(maximum_log, term_log)
      end do
      value = exp(maximum_log) * sum(exp(logs - maximum_log))
      value = min(max(value, 0.0_dp), 1.0_dp)
   end function binomial_cdf

end module quarks_stats
