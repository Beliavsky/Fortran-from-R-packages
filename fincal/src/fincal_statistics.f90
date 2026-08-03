! SPDX-License-Identifier: GPL-2.0-or-later
module fincal_statistics
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use fincal_kinds, only : dp
   use fincal_status, only : fincal_ok, fincal_invalid_input, fincal_size_mismatch, fincal_weights_not_unit
   use fincal_rates, only : hpr
   implicit none
   private

   public :: geometric_mean, harmonic_mean, sampling_error
   public :: weighted_portfolio_return, twrr
   public :: wpr

   interface wpr
      module procedure weighted_portfolio_return
   end interface wpr
contains
   function geometric_mean(returns, status) result(value)
      real(dp), intent(in) :: returns(:)
      integer, intent(out), optional :: status
      real(dp) :: value

      if (size(returns) == 0 .or. any(returns <= -1.0_dp)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincal_invalid_input
         return
      end if
      value = exp(sum(log(1.0_dp + returns)) / real(size(returns), dp)) - 1.0_dp
      if (present(status)) status = fincal_ok
   end function geometric_mean

   function harmonic_mean(values, status) result(value)
      real(dp), intent(in) :: values(:)
      integer, intent(out), optional :: status
      real(dp) :: value

      if (size(values) == 0 .or. any(abs(values) <= tiny(1.0_dp))) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincal_invalid_input
         return
      end if
      value = real(size(values), dp) / sum(1.0_dp / values)
      if (present(status)) status = fincal_ok
   end function harmonic_mean

   elemental pure function sampling_error(sample_mean, population_mean) result(value)
      real(dp), intent(in) :: sample_mean, population_mean
      real(dp) :: value
      value = sample_mean - population_mean
   end function sampling_error

   function weighted_portfolio_return(returns, weights, status) result(value)
      real(dp), intent(in) :: returns(:), weights(:)
      integer, intent(out), optional :: status
      real(dp) :: value

      if (size(returns) /= size(weights)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincal_size_mismatch
         return
      end if
      value = sum(returns * weights)
      if (present(status)) then
         if (abs(sum(weights) - 1.0_dp) <= 100.0_dp * epsilon(1.0_dp)) then
            status = fincal_ok
         else
            status = fincal_weights_not_unit
         end if
      end if
   end function weighted_portfolio_return

   function twrr(ending_values, beginning_values, cash_flows_received, status) result(value)
      real(dp), intent(in) :: ending_values(:), beginning_values(:), cash_flows_received(:)
      integer, intent(out), optional :: status
      real(dp) :: value
      integer :: i, n

      n = size(ending_values)
      if (n == 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincal_invalid_input
         return
      end if
      if (size(beginning_values) /= n .or. size(cash_flows_received) /= n) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincal_size_mismatch
         return
      end if

      value = 1.0_dp
      do i = 1, n
         value = value * (1.0_dp + hpr(ending_values(i), beginning_values(i), cash_flows_received(i)))
      end do
      value = value ** (1.0_dp / real(n, dp)) - 1.0_dp
      if (present(status)) status = fincal_ok
   end function twrr
end module fincal_statistics
