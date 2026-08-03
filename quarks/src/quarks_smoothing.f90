module quarks_smoothing
   use quarks_kinds, only : dp
   use quarks_types, only : smooth_none, smooth_lpr, smooth_auto
   implicit none
   private
   public :: smooth_scale

contains

   subroutine smooth_scale(x, mode, x_standardized, forecast_scale, bandwidth, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: mode
      real(dp), allocatable, intent(out) :: x_standardized(:)
      real(dp), intent(out) :: forecast_scale
      real(dp), intent(in), optional :: bandwidth
      integer, intent(out), optional :: status
      real(dp), allocatable :: local_variance(:)
      real(dp) :: b, u, v, weight, numerator, denominator, global_variance
      integer :: i, j, n
      n = size(x)
      allocate(x_standardized(n), local_variance(n))
      if (present(status)) status = 0
      if (n <= 1) then
         x_standardized = x
         forecast_scale = 1.0_dp
         if (present(status)) status = 1
         return
      end if
      if (mode == smooth_none) then
         x_standardized = x
         forecast_scale = 1.0_dp
         return
      end if
      if (present(bandwidth)) then
         b = bandwidth
      else if (mode == smooth_lpr) then
         b = 0.15_dp
      else if (mode == smooth_auto) then
         b = min(0.30_dp, max(0.06_dp, 1.06_dp * real(n, dp)**(-0.2_dp)))
      else
         x_standardized = x
         forecast_scale = 1.0_dp
         if (present(status)) status = 1
         return
      end if
      b = max(b, 1.0_dp / real(max(2, n - 1), dp))
      global_variance = sum(x * x) / real(n, dp)
      global_variance = max(global_variance, 1.0e-16_dp)
      do i = 1, n
         u = real(i - 1, dp) / real(n - 1, dp)
         numerator = 0.0_dp
         denominator = 0.0_dp
         do j = 1, n
            v = real(j - 1, dp) / real(n - 1, dp)
            weight = exp(-0.5_dp * ((u - v) / b)**2)
            numerator = numerator + weight * x(j)**2
            denominator = denominator + weight
         end do
         local_variance(i) = numerator / max(denominator, tiny(1.0_dp))
         local_variance(i) = max(local_variance(i), 1.0e-12_dp * global_variance)
      end do
      x_standardized = x / sqrt(local_variance)
      u = 1.0_dp + 1.0_dp / real(n - 1, dp)
      numerator = 0.0_dp
      denominator = 0.0_dp
      do j = 1, n
         v = real(j - 1, dp) / real(n - 1, dp)
         weight = exp(-0.5_dp * ((u - v) / b)**2)
         numerator = numerator + weight * x(j)**2
         denominator = denominator + weight
      end do
      forecast_scale = sqrt(max(numerator / max(denominator, tiny(1.0_dp)), &
         1.0e-12_dp * global_variance))
   end subroutine smooth_scale

end module quarks_smoothing
