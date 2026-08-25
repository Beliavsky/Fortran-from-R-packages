! SPDX-License-Identifier: GPL-2.0-or-later
module fints_summary_mod
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use fints_kinds, only : dp
   use fints_status, only : fints_ok, fints_no_data, fints_invalid_input
   use fints_types, only : summary_result
   use r_descriptive, only : r_count_nonmissing, r_mean, r_variance
   implicit none
   private
   public :: fints_summary, sample_mean, sample_variance, finite_count

contains

   pure integer function finite_count(x) result(n)
      real(dp), intent(in) :: x(:)
      n = r_count_nonmissing(x)
   end function finite_count

   pure real(dp) function sample_mean(x) result(mean_value)
      real(dp), intent(in) :: x(:)
      mean_value = r_mean(x, na_rm=.true.)
      if (ieee_is_nan(mean_value)) mean_value = 0.0_dp
   end function sample_mean

   pure real(dp) function sample_variance(x, mean_value) result(variance)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: mean_value
      real(dp) :: center

      if (present(mean_value)) then
         center = mean_value
      else
         center = sample_mean(x)
      end if
      variance = r_variance(x, na_rm=.true., center=center)
      if (ieee_is_nan(variance)) variance = 0.0_dp
   end function sample_variance

   subroutine fints_summary(x, result, start)
      real(dp), intent(in) :: x(:)
      type(summary_result), intent(out) :: result
      real(dp), intent(in), optional :: start
      real(dp) :: mean_value, m2, m3, m4, dx
      integer :: i, n

      result = summary_result()
      if (present(start)) result%start = start
      n = finite_count(x)
      result%size = n
      if (n < 1) then
         result%status = fints_no_data
         return
      end if

      mean_value = sample_mean(x)
      result%mean = mean_value
      result%minimum = huge(1.0_dp)
      result%maximum = -huge(1.0_dp)
      m2 = 0.0_dp
      m3 = 0.0_dp
      m4 = 0.0_dp
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) cycle
         dx = x(i) - mean_value
         m2 = m2 + dx ** 2
         m3 = m3 + dx ** 3
         m4 = m4 + dx ** 4
         result%minimum = min(result%minimum, x(i))
         result%maximum = max(result%maximum, x(i))
      end do

      if (n > 1) result%standard_deviation = sqrt(m2 / real(n - 1, dp))
      if (m2 > 0.0_dp) then
         result%skewness = (m3 / real(n, dp)) / (m2 / real(n, dp)) ** 1.5_dp
         result%excess_kurtosis = (m4 / real(n, dp)) / (m2 / real(n, dp)) ** 2 - 3.0_dp
      else
         result%skewness = 0.0_dp
         result%excess_kurtosis = 0.0_dp
      end if
      result%status = fints_ok
   end subroutine fints_summary

end module fints_summary_mod
