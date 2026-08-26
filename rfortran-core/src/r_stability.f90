! SPDX-License-Identifier: MIT
module r_stability
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_negative_inf, ieee_positive_inf, &
      ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   implicit none
   private

   public :: r_log_mean_exp, r_log_sum_exp

contains

   pure real(dp) function r_log_sum_exp(x, na_rm) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: na_rm
      real(dp) :: maximum
      logical :: remove_na
      integer :: i, n

      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      if (.not. remove_na .and. any(ieee_is_nan(x))) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      n = count(.not. ieee_is_nan(x))
      if (n == 0) then
         value = ieee_value(0.0_dp, ieee_negative_inf)
         return
      end if
      maximum = -huge(1.0_dp)
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) maximum = max(maximum, x(i))
      end do
      if (maximum == ieee_value(0.0_dp, ieee_positive_inf)) then
         value = maximum
      else if (maximum == ieee_value(0.0_dp, ieee_negative_inf)) then
         value = maximum
      else
         value = maximum + log(sum(exp(x - maximum), mask=.not. ieee_is_nan(x)))
      end if
   end function r_log_sum_exp

   pure real(dp) function r_log_mean_exp(x, na_rm) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: na_rm
      logical :: remove_na
      integer :: n

      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      if (.not. remove_na .and. any(ieee_is_nan(x))) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      n = count(.not. ieee_is_nan(x))
      if (n == 0) then
         value = ieee_value(0.0_dp, ieee_negative_inf)
      else
         value = r_log_sum_exp(x, remove_na) - log(real(n, dp))
      end if
   end function r_log_mean_exp

end module r_stability
