! SPDX-License-Identifier: MIT
! SPDX-FileComment: Classical seasonal decomposition compatible with R stats conventions.
module r_stats_decomposition
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
   use r_kinds, only: dp
   use r_stats_types, only: seasonal_decomposition_t
   implicit none
   private

   public :: decompose

contains

   pure function decompose(x, period, multiplicative, filter) result(decomposition)
      !! Decomposes a complete finite series using R's default centered moving average.
      real(dp), intent(in) :: x(:) !! Finite observations spanning at least two cycles.
      integer, intent(in) :: period !! Number of observations in one seasonal cycle.
      logical, intent(in), optional :: multiplicative !! Select multiplicative decomposition.
      real(dp), intent(in), optional :: filter(:) !! Centered convolution weights.
      type(seasonal_decomposition_t) :: decomposition
      real(dp) :: mean_figure, total
      real(dp), allocatable :: preliminary(:), weights(:)
      integer :: filter_length, first_valid, i, j, last_valid, observations
      integer :: seasonal_position, shift
      logical :: use_multiplicative

      use_multiplicative = .false.
      if (present(multiplicative)) use_multiplicative = multiplicative
      decomposition%period = period
      decomposition%multiplicative = use_multiplicative
      if (period <= 1 .or. size(x) < 2*period) then
         decomposition%status = 1
         return
      end if
      if (.not. all(ieee_is_finite(x))) then
         decomposition%status = 2
         return
      end if

      if (present(filter)) then
         if (size(filter) < 1 .or. size(filter) > size(x) .or. &
             .not. all(ieee_is_finite(filter))) then
            decomposition%status = 2
            return
         end if
         weights = filter
      else if (modulo(period, 2) == 0) then
         allocate (weights(period + 1), source=1.0_dp/real(period, dp))
         weights([1, period + 1]) = 0.5_dp/real(period, dp)
      else
         allocate (weights(period), source=1.0_dp/real(period, dp))
      end if

      filter_length = size(weights)
      shift = filter_length/2
      first_valid = filter_length - shift
      last_valid = size(x) - shift
      allocate (decomposition%trend(size(x)), decomposition%seasonal(size(x)), &
                decomposition%random(size(x)), decomposition%figure(period), &
                preliminary(size(x)))
      decomposition%trend = ieee_value(0.0_dp, ieee_quiet_nan)
      preliminary = ieee_value(0.0_dp, ieee_quiet_nan)
      decomposition%random = ieee_value(0.0_dp, ieee_quiet_nan)

      do i = first_valid, last_valid
         total = 0.0_dp
         do j = 1, filter_length
            total = total + weights(j)*x(i + shift - j + 1)
         end do
         decomposition%trend(i) = total
         if (use_multiplicative) then
            if (total == 0.0_dp) then
               decomposition%status = 2
               return
            end if
            preliminary(i) = x(i)/total
         else
            preliminary(i) = x(i) - total
         end if
      end do

      do seasonal_position = 1, period
         total = 0.0_dp
         observations = 0
         do i = first_valid, last_valid
            if (modulo(i - 1, period) + 1 /= seasonal_position) cycle
            total = total + preliminary(i)
            observations = observations + 1
         end do
         if (observations == 0) then
            decomposition%status = 3
            return
         end if
         decomposition%figure(seasonal_position) = total/real(observations, dp)
      end do
      mean_figure = sum(decomposition%figure)/real(period, dp)
      if (use_multiplicative) then
         if (mean_figure == 0.0_dp) then
            decomposition%status = 2
            return
         end if
         decomposition%figure = decomposition%figure/mean_figure
      else
         decomposition%figure = decomposition%figure - mean_figure
      end if

      do i = 1, size(x)
         seasonal_position = modulo(i - 1, period) + 1
         decomposition%seasonal(i) = decomposition%figure(seasonal_position)
      end do
      do i = first_valid, last_valid
         if (use_multiplicative) then
            decomposition%random(i) = x(i)/decomposition%seasonal(i)/decomposition%trend(i)
         else
            decomposition%random(i) = x(i) - decomposition%seasonal(i) - decomposition%trend(i)
         end if
      end do
   end function decompose

end module r_stats_decomposition
