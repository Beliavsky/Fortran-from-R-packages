! SPDX-License-Identifier: GPL-2.0-or-later
module fints_finance
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use fints_kinds, only : dp
   implicit none
   private
   public :: compound_interest, simple_to_log_returns

   interface compound_interest
      module procedure compound_interest_scalar
      module procedure compound_interest_vector
   end interface compound_interest

   interface simple_to_log_returns
      module procedure simple_to_log_returns_scalar
      module procedure simple_to_log_returns_vector
   end interface simple_to_log_returns

contains

   pure real(dp) function compound_interest_scalar(interest, periods, frequency, &
      net_value) result(value)
      real(dp), intent(in) :: interest
      real(dp), intent(in), optional :: periods, frequency
      logical, intent(in), optional :: net_value
      real(dp) :: number_periods, compounding
      logical :: return_net_increase

      number_periods = 1.0_dp
      if (present(periods)) number_periods = periods
      compounding = 1.0_dp
      if (present(frequency)) compounding = frequency
      if (ieee_is_finite(compounding)) then
         value = exp(number_periods * compounding * log_one_plus(interest / compounding))
      else
         value = exp(number_periods * interest)
      end if
      return_net_increase = .false.
      if (present(net_value)) return_net_increase = net_value
      if (return_net_increase) value = value - 1.0_dp
   end function compound_interest_scalar

   pure function compound_interest_vector(interest, periods, frequency, &
      net_value) result(value)
      real(dp), intent(in) :: interest(:)
      real(dp), intent(in), optional :: periods(:), frequency(:)
      logical, intent(in), optional :: net_value
      real(dp) :: value(size(interest))
      real(dp) :: p, f
      integer :: i

      do i = 1, size(interest)
         p = 1.0_dp
         f = 1.0_dp
         if (present(periods)) p = periods(min(i, size(periods)))
         if (present(frequency)) f = frequency(min(i, size(frequency)))
         value(i) = compound_interest_scalar(interest(i), p, f, net_value)
      end do
   end function compound_interest_vector

   pure real(dp) function simple_to_log_returns_scalar(simple_return) result(log_return)
      real(dp), intent(in) :: simple_return
      log_return = log_one_plus(simple_return)
   end function simple_to_log_returns_scalar

   pure function simple_to_log_returns_vector(simple_return) result(log_return)
      real(dp), intent(in) :: simple_return(:)
      real(dp) :: log_return(size(simple_return))
      log_return = log_one_plus(simple_return)
   end function simple_to_log_returns_vector

   pure elemental real(dp) function log_one_plus(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: term
      integer :: k

      if (abs(x) > 1.0e-6_dp) then
         value = log(1.0_dp + x)
      else
         value = 0.0_dp
         term = x
         do k = 1, 12
            value = value + term / real(k, dp)
            term = -term * x
         end do
      end if
   end function log_one_plus

end module fints_finance
