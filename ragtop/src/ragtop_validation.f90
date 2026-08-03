! SPDX-License-Identifier: GPL-2.0-or-later
module ragtop_validation
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ragtop_kinds, only : dp
   use ragtop_constants, only : ragtop_ok, ragtop_invalid_argument
   use ragtop_types, only : market_spec
   use ragtop_term_structures, only : discount_factor, cumulative_variance, &
                                      survival_probability, default_intensity
   implicit none
   private
   public :: check_discount_factor, check_variance_cumulation
   public :: check_survival_probability, check_default_intensity

contains

   subroutine check_discount_factor(market, time_points, status)
      type(market_spec), intent(in) :: market
      real(dp), intent(in) :: time_points(:)
      integer, intent(out) :: status
      real(dp) :: previous, current
      integer :: i
      status = ragtop_ok
      if (any(time_points < 0.0_dp)) then
         status = ragtop_invalid_argument
         return
      end if
      previous = 1.0_dp
      do i = 1, size(time_points)
         current = discount_factor(market,time_points(i),0.0_dp)
         if (.not. ieee_is_finite(current) .or. current <= 0.0_dp .or. &
             current > previous+100.0_dp*epsilon(1.0_dp)) then
            status = ragtop_invalid_argument
            return
         end if
         previous = current
      end do
   end subroutine check_discount_factor

   subroutine check_variance_cumulation(market, time_points, status)
      type(market_spec), intent(in) :: market
      real(dp), intent(in) :: time_points(:)
      integer, intent(out) :: status
      real(dp) :: previous, current
      integer :: i
      status = ragtop_ok
      previous = 0.0_dp
      do i = 1, size(time_points)
         current = cumulative_variance(market,time_points(i),0.0_dp)
         if (.not. ieee_is_finite(current) .or. current < previous .or. &
             current < 0.0_dp) then
            status = ragtop_invalid_argument
            return
         end if
         previous = current
      end do
   end subroutine check_variance_cumulation

   subroutine check_survival_probability(market, time_points, status)
      type(market_spec), intent(in) :: market
      real(dp), intent(in) :: time_points(:)
      integer, intent(out) :: status
      real(dp) :: previous, current
      integer :: i
      status = ragtop_ok
      previous = 1.0_dp
      do i = 1, size(time_points)
         current = survival_probability(market,time_points(i),0.0_dp)
         if (.not. ieee_is_finite(current) .or. current < 0.0_dp .or. &
             current > 1.0_dp .or. &
             current > previous+100.0_dp*epsilon(1.0_dp)) then
            status = ragtop_invalid_argument
            return
         end if
         previous = current
      end do
   end subroutine check_survival_probability

   subroutine check_default_intensity(market, time, stock, status)
      type(market_spec), intent(in) :: market
      real(dp), intent(in) :: time, stock(:)
      integer, intent(out) :: status
      real(dp), allocatable :: intensity(:)
      allocate(intensity(size(stock)))
      intensity = default_intensity(market,time,stock)
      if (all(ieee_is_finite(intensity)) .and. all(intensity >= 0.0_dp)) then
         status = ragtop_ok
      else
         status = ragtop_invalid_argument
      end if
   end subroutine check_default_intensity

end module ragtop_validation
