! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Dmitriy Mayorov
module vasicekfit_normal
   use, intrinsic :: iso_fortran_env, only : int64
   use vasicekfit_kinds, only : dp, pi, sqrt_two, log_two_pi
   implicit none
   private

   public :: normal_pdf, normal_cdf, normal_quantile
   public :: fill_standard_normals, seed_random_number

contains

   pure elemental real(dp) function normal_pdf(x) result(value)
      real(dp), intent(in) :: x
      value = exp(-0.5_dp * x * x - 0.5_dp * log_two_pi)
   end function normal_pdf

   pure elemental real(dp) function normal_cdf(x) result(value)
      real(dp), intent(in) :: x
      value = 0.5_dp * erfc(-x / sqrt_two)
   end function normal_cdf

   pure elemental real(dp) function normal_quantile(probability) result(value)
      real(dp), intent(in) :: probability
      real(dp) :: q, r
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e+01_dp,  2.209460984245205e+02_dp, &
         -2.759285104469687e+02_dp,  1.383577518672690e+02_dp, &
         -3.066479806614716e+01_dp,  2.506628277459239e+00_dp ]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e+01_dp,  1.615858368580409e+02_dp, &
         -1.556989798598866e+02_dp,  6.680131188771972e+01_dp, &
         -1.328068155288572e+01_dp ]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-03_dp, -3.223964580411365e-01_dp, &
         -2.400758277161838e+00_dp, -2.549732539343734e+00_dp, &
          4.374664141464968e+00_dp,  2.938163982698783e+00_dp ]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462e-03_dp,  3.224671290700398e-01_dp, &
          2.445134137142996e+00_dp,  3.754408661907416e+00_dp ]
      real(dp), parameter :: p_low = 0.02425_dp
      real(dp), parameter :: p_high = 1.0_dp - p_low

      if (probability <= 0.0_dp) then
         value = -huge(1.0_dp)
      else if (probability >= 1.0_dp) then
         value = huge(1.0_dp)
      else if (probability < p_low) then
         q = sqrt(-2.0_dp * log(probability))
         value = (((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + c(5)) * q + c(6)) / &
            ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0_dp)
      else if (probability <= p_high) then
         q = probability - 0.5_dp
         r = q * q
         value = (((((a(1) * r + a(2)) * r + a(3)) * r + a(4)) * r + a(5)) * r + a(6)) * q / &
            (((((b(1) * r + b(2)) * r + b(3)) * r + b(4)) * r + b(5)) * r + 1.0_dp)
      else
         q = sqrt(-2.0_dp * log(1.0_dp - probability))
         value = -(((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + c(5)) * q + c(6)) / &
            ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0_dp)
      end if

      if (probability > 0.0_dp .and. probability < 1.0_dp) then
         value = value - (normal_cdf(value) - probability) / normal_pdf(value)
      end if
   end function normal_quantile

   subroutine seed_random_number(seed)
      integer(int64), intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)
      integer(int64) :: state

      call random_seed(size=n)
      allocate(put(n))
      state = seed
      do i = 1, n
         state = modulo(6364136223846793005_int64 * state + 1442695040888963407_int64, &
            huge(1_int64))
         put(i) = int(modulo(state, int(huge(1), int64)), kind(put(i)))
         if (put(i) == 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine seed_random_number

   subroutine fill_standard_normals(values, seed)
      real(dp), intent(out) :: values(:)
      integer(int64), intent(in), optional :: seed
      integer :: i
      real(dp) :: u1, u2, radius, angle

      if (present(seed)) call seed_random_number(seed)
      i = 1
      do while (i <= size(values))
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1, tiny(1.0_dp))
         radius = sqrt(-2.0_dp * log(u1))
         angle = 2.0_dp * pi * u2
         values(i) = radius * cos(angle)
         if (i + 1 <= size(values)) values(i + 1) = radius * sin(angle)
         i = i + 2
      end do
   end subroutine fill_standard_normals

end module vasicekfit_normal
