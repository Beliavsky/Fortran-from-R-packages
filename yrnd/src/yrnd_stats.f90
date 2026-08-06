! SPDX-License-Identifier: GPL-3.0-only
module yrnd_stats
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use yrnd_kinds, only : dp
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)

   public :: normal_cdf, normal_quantile, lognormal_pdf, lognormal_cdf
   public :: lognormal_quantile, trapezoid_integral, normalize_density
   public :: seed_random, random_normal, sample_mean, sample_sd, gaussian_kde

contains

   elemental real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function normal_cdf

   elemental real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
         -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
         -3.066479806614716e1_dp, 2.506628277459239_dp ]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
         -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
         -1.328068155288572e1_dp ]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
         -2.400758277161838_dp, -2.549732539343734_dp, &
         4.374664141464968_dp, 2.938163982698783_dp ]
      real(dp), parameter :: d(4) = [ &
         7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
         2.445134137142996_dp, 3.754408661907416_dp ]
      real(dp), parameter :: plow = 0.02425_dp
      real(dp), parameter :: phigh = 1.0_dp - plow
      real(dp) :: q, r

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp * log(p))
         x = (((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + c(5)) * q + c(6)) / &
            ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0_dp)
      else if (p <= phigh) then
         q = p - 0.5_dp
         r = q * q
         x = (((((a(1) * r + a(2)) * r + a(3)) * r + a(4)) * r + a(5)) * r + a(6)) * q / &
            (((((b(1) * r + b(2)) * r + b(3)) * r + b(4)) * r + b(5)) * r + 1.0_dp)
      else
         q = sqrt(-2.0_dp * log(1.0_dp - p))
         x = -(((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + c(5)) * q + c(6)) / &
            ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0_dp)
      end if
   end function normal_quantile

   elemental real(dp) function lognormal_pdf(x, meanlog, sdlog) result(f)
      real(dp), intent(in) :: x, meanlog, sdlog
      if (x <= 0.0_dp .or. sdlog <= 0.0_dp) then
         f = 0.0_dp
      else
         f = exp(-0.5_dp * ((log(x) - meanlog) / sdlog) ** 2) / &
            (x * sdlog * sqrt(2.0_dp * pi))
      end if
   end function lognormal_pdf

   elemental real(dp) function lognormal_cdf(x, meanlog, sdlog) result(p)
      real(dp), intent(in) :: x, meanlog, sdlog
      if (x <= 0.0_dp) then
         p = 0.0_dp
      else if (sdlog <= 0.0_dp) then
         p = merge(1.0_dp, 0.0_dp, x >= exp(meanlog))
      else
         p = normal_cdf((log(x) - meanlog) / sdlog)
      end if
   end function lognormal_cdf

   elemental real(dp) function lognormal_quantile(p, meanlog, sdlog) result(x)
      real(dp), intent(in) :: p, meanlog, sdlog
      if (p <= 0.0_dp) then
         x = 0.0_dp
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else
         x = exp(meanlog + sdlog * normal_quantile(p))
      end if
   end function lognormal_quantile

   pure real(dp) function trapezoid_integral(x, y) result(value)
      real(dp), intent(in) :: x(:), y(:)
      integer :: i
      if (size(x) /= size(y) .or. size(x) < 2) then
         value = 0.0_dp
         return
      end if
      value = 0.0_dp
      do i = 1, size(x) - 1
         value = value + 0.5_dp * (y(i) + y(i + 1)) * (x(i + 1) - x(i))
      end do
   end function trapezoid_integral

   subroutine normalize_density(x, density)
      real(dp), intent(in) :: x(:)
      real(dp), intent(inout) :: density(:)
      real(dp) :: area
      area = trapezoid_integral(x, density)
      if (area > 0.0_dp) density = density / area
   end subroutine normalize_density

   subroutine seed_random(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: values(:)
      call random_seed(size=n)
      allocate(values(n))
      do i = 1, n
         values(i) = modulo(seed + 104729 * i, huge(1) - 1)
         if (values(i) <= 0) values(i) = i
      end do
      call random_seed(put=values)
   end subroutine seed_random

   real(dp) function random_normal() result(z)
      real(dp) :: u1, u2
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
   end function random_normal

   pure real(dp) function sample_mean(x) result(value)
      real(dp), intent(in) :: x(:)
      if (size(x) == 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         value = sum(x) / real(size(x), dp)
      end if
   end function sample_mean

   pure real(dp) function sample_sd(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      if (size(x) < 2) then
         value = 0.0_dp
      else
         m = sample_mean(x)
         value = sqrt(sum((x - m) ** 2) / real(size(x) - 1, dp))
      end if
   end function sample_sd

   subroutine gaussian_kde(samples, grid, density, bandwidth)
      real(dp), intent(in) :: samples(:), grid(:)
      real(dp), intent(out) :: density(size(grid))
      real(dp), intent(in), optional :: bandwidth
      real(dp) :: h, scale
      integer :: i
      if (size(samples) == 0) then
         density = 0.0_dp
         return
      end if
      if (present(bandwidth)) then
         h = bandwidth
      else
         h = 1.06_dp * max(sample_sd(samples), sqrt(epsilon(1.0_dp))) * &
            real(size(samples), dp) ** (-0.2_dp)
      end if
      h = max(h, sqrt(epsilon(1.0_dp)))
      scale = 1.0_dp / (real(size(samples), dp) * h * sqrt(2.0_dp * pi))
      do i = 1, size(grid)
         density(i) = scale * sum(exp(-0.5_dp * ((grid(i) - samples) / h) ** 2))
      end do
      call normalize_density(grid, density)
   end subroutine gaussian_kde

end module yrnd_stats
