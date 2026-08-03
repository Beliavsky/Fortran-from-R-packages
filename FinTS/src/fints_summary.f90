! SPDX-License-Identifier: GPL-2.0-or-later
module fints_summary_mod
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use fints_kinds, only : dp
   use fints_status, only : fints_ok, fints_no_data, fints_invalid_input
   use fints_types, only : summary_result
   implicit none
   private
   public :: fints_summary, sample_mean, sample_variance, finite_count

contains

   pure integer function finite_count(x) result(n)
      real(dp), intent(in) :: x(:)
      integer :: i

      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) n = n + 1
      end do
   end function finite_count

   pure real(dp) function sample_mean(x) result(mean_value)
      real(dp), intent(in) :: x(:)
      real(dp) :: total
      integer :: i, n

      total = 0.0_dp
      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) then
            total = total + x(i)
            n = n + 1
         end if
      end do
      if (n > 0) then
         mean_value = total / real(n, dp)
      else
         mean_value = 0.0_dp
      end if
   end function sample_mean

   pure real(dp) function sample_variance(x, mean_value) result(variance)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: mean_value
      real(dp) :: center, total
      integer :: i, n

      if (present(mean_value)) then
         center = mean_value
      else
         center = sample_mean(x)
      end if
      total = 0.0_dp
      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) then
            total = total + (x(i) - center) ** 2
            n = n + 1
         end if
      end do
      if (n > 1) then
         variance = total / real(n - 1, dp)
      else
         variance = 0.0_dp
      end if
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
