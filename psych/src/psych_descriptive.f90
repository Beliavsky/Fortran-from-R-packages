module psych_descriptive
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use psych_stats, only : finite_count, finite_mean, finite_sd, finite_median
   use psych_basic, only : skewness, excess_kurtosis
   use psych_types, only : describe_result
   implicit none
   private
   public :: describe_vector, describe_matrix, top_bottom_means, response_frequencies
   public :: pairwise_count_matrix, pairwise_zero_matrix, pairwise_impute_mean

contains

   pure subroutine describe_vector(x, result)
      real(dp), intent(in) :: x(:) !! Numeric vector summarized using psych-style descriptive statistics.
      type(describe_result), intent(out) :: result !! Count, center, spread, range, skewness, kurtosis, and standard error.
      integer :: i
      result%n = finite_count(x)
      result%mean = finite_mean(x)
      result%sd = finite_sd(x)
      result%median = finite_median(x)
      result%minimum = huge(1.0_dp)
      result%maximum = -huge(1.0_dp)
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) then
            result%minimum = min(result%minimum,x(i))
            result%maximum = max(result%maximum,x(i))
         end if
      end do
      if (result%n == 0) then
         result%minimum = 0.0_dp
         result%maximum = 0.0_dp
      end if
      result%range = result%maximum-result%minimum
      result%skew = skewness(x,3)
      result%kurtosis = excess_kurtosis(x,3)
      if (result%n > 0) then
         result%se = result%sd/sqrt(real(result%n,dp))
      else
         result%se = 0.0_dp
      end if
   end subroutine describe_vector

   pure subroutine describe_matrix(x, results)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix summarized columnwise.
      type(describe_result), intent(out) :: results(:) !! One descriptive-statistics result per input column.
      integer :: j
      do j = 1, min(size(x,2),size(results))
         call describe_vector(x(:,j),results(j))
      end do
   end subroutine describe_matrix

   pure subroutine top_bottom_means(x, fraction, low_mean, high_mean)
      real(dp), intent(in) :: x(:) !! Numeric observations split by lower and upper quantile cutoffs.
      real(dp), intent(in), optional :: fraction !! Fraction in each tail; default `0.27`.
      real(dp), intent(out) :: low_mean !! Mean of observations in the lower tail.
      real(dp), intent(out) :: high_mean !! Mean of observations in the upper tail.
      real(dp), allocatable :: work(:)
      real(dp) :: frac
      integer :: i, k, n
      frac = 0.27_dp
      if (present(fraction)) frac = max(0.0_dp,min(0.5_dp,fraction))
      n = finite_count(x)
      if (n == 0) then
         low_mean = 0.0_dp
         high_mean = 0.0_dp
         return
      end if
      allocate(work(n))
      k = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) then
            k = k+1
            work(k) = x(i)
         end if
      end do
      call sort_local(work)
      k = max(1,nint(frac*real(n,dp)))
      low_mean = sum(work(1:k))/real(k,dp)
      high_mean = sum(work(n-k+1:n))/real(k,dp)
   end subroutine top_bottom_means

   pure subroutine response_frequencies(x, minimum, maximum, counts)
      integer, intent(in) :: x(:) !! Integer-coded item responses.
      integer, intent(in) :: minimum !! Lowest response category represented in `counts(1)`.
      integer, intent(in) :: maximum !! Highest response category represented in the output.
      integer, intent(out) :: counts(:) !! Category counts for values from `minimum` through `maximum`.
      integer :: i, index
      counts = 0
      do i = 1, size(x)
         if (x(i) >= minimum .and. x(i) <= maximum) then
            index = x(i)-minimum+1
            if (index <= size(counts)) counts(index) = counts(index)+1
         end if
      end do
   end subroutine response_frequencies

   pure subroutine pairwise_count_matrix(x, counts)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix whose NaNs represent missing data.
      integer, intent(out) :: counts(:, :) !! Number of rows observed for every variable pair.
      integer :: i, j, k, n
      do j = 1, size(x,2)
         do i = 1, size(x,2)
            n = 0
            do k = 1, size(x,1)
               if (.not. ieee_is_nan(x(k,i)) .and. .not. ieee_is_nan(x(k,j))) n = n+1
            end do
            counts(i,j) = n
         end do
      end do
   end subroutine pairwise_count_matrix

   pure subroutine pairwise_zero_matrix(x, result)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix whose NaNs represent missing data.
      real(dp), intent(out) :: result(:, :) !! Copy of `x` with missing values replaced by zero.
      integer :: i, j
      result = x
      do j = 1, size(x,2)
         do i = 1, size(x,1)
            if (ieee_is_nan(result(i,j))) result(i,j) = 0.0_dp
         end do
      end do
   end subroutine pairwise_zero_matrix

   pure subroutine pairwise_impute_mean(x, result)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix whose NaNs represent missing data.
      real(dp), intent(out) :: result(:, :) !! Matrix with each NaN replaced by its column mean.
      real(dp) :: m
      integer :: i, j
      result = x
      do j = 1, size(x,2)
         m = finite_mean(x(:,j))
         do i = 1, size(x,1)
            if (ieee_is_nan(result(i,j))) result(i,j) = m
         end do
      end do
   end subroutine pairwise_impute_mean

   pure subroutine sort_local(x)
      real(dp), intent(inout) :: x(:) !! Vector sorted into ascending order in place.
      real(dp) :: key
      integer :: i, j
      do i = 2, size(x)
         key = x(i)
         j = i-1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j+1) = x(j)
            j = j-1
         end do
         x(j+1) = key
      end do
   end subroutine sort_local

end module psych_descriptive
