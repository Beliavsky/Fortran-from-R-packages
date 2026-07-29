! SPDX-License-Identifier: GPL-2.0-only
!
! Modern Fortran translation of computational methods from GCPM 1.2.2.
! Original software copyright (C) 2015 Kevin Jakob and Dr. Matthias Fischer.
! Fortran translation copyright (C) 2026.
module gcpm_math
   use, intrinsic :: iso_fortran_env, only: int64
   use gcpm_kinds, only: dp
   implicit none
   private

   public :: normal_cdf
   public :: normal_quantile
   public :: seed_random_number
   public :: random_poisson
   public :: random_standard_normal
   public :: random_gamma
   public :: correlation_matrix
   public :: quadratic_form

contains

   pure elemental function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      real(dp) :: p

      p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function normal_cdf

   pure elemental function normal_quantile(p) result(x)
      ! Peter J. Acklam's rational approximation, followed by one
      ! Halley refinement. Accurate to near double precision.
      real(dp), intent(in) :: p
      real(dp) :: x
      real(dp) :: q, r, e, u
      real(dp), parameter :: p_low = 0.02425_dp
      real(dp), parameter :: p_high = 1.0_dp - p_low
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376d+01, 2.209460984245205d+02, &
         -2.759285104469687d+02, 1.383577518672690d+02, &
         -3.066479806614716d+01, 2.506628277459239d+00 ]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406d+01, 1.615858368580409d+02, &
         -1.556989798598866d+02, 6.680131188771972d+01, &
         -1.328068155288572d+01 ]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293d-03, -3.223964580411365d-01, &
         -2.400758277161838d+00, -2.549732539343734d+00, &
          4.374664141464968d+00,  2.938163982698783d+00 ]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462d-03, 3.224671290700398d-01, &
          2.445134137142996d+00, 3.754408661907416d+00 ]

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
         return
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
         return
      end if

      if (p < p_low) then
         q = sqrt(-2.0_dp * log(p))
         x = (((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + c(5)) * q + c(6)) / &
             ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0_dp)
      else if (p <= p_high) then
         q = p - 0.5_dp
         r = q * q
         x = (((((a(1) * r + a(2)) * r + a(3)) * r + a(4)) * r + a(5)) * r + a(6)) * q / &
             (((((b(1) * r + b(2)) * r + b(3)) * r + b(4)) * r + b(5)) * r + 1.0_dp)
      else
         q = sqrt(-2.0_dp * log(1.0_dp - p))
         x = -(((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + c(5)) * q + c(6)) / &
              ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0_dp)
      end if

      e = normal_cdf(x) - p
      u = e * sqrt(2.0_dp * acos(-1.0_dp)) * exp(0.5_dp * x * x)
      x = x - u / (1.0_dp + 0.5_dp * x * u)
   end function normal_quantile

   subroutine seed_random_number(seed)
      integer, intent(in) :: seed
      integer :: i, n
      integer, allocatable :: put(:)
      integer(int64) :: value

      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         value = int(abs(seed), int64) + 104729_int64 * int(i, int64) + &
                 12345_int64 * int(i * i, int64)
         put(i) = int(mod(value, 2147483646_int64)) + 1
      end do
      call random_seed(put=put)
   end subroutine seed_random_number

   function random_standard_normal() result(z)
      real(dp) :: z
      real(dp) :: u1, u2

      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)
   end function random_standard_normal

   recursive function random_gamma(shape, scale) result(x)
      ! Marsaglia-Tsang gamma generator.
      real(dp), intent(in) :: shape
      real(dp), intent(in), optional :: scale
      real(dp) :: x
      real(dp) :: a, d, c, u, v, z, scl

      scl = 1.0_dp
      if (present(scale)) scl = scale
      if (shape <= 0.0_dp .or. scl <= 0.0_dp) then
         x = 0.0_dp
         return
      end if

      if (shape < 1.0_dp) then
         call random_number(u)
         x = random_gamma(shape + 1.0_dp, scl) * u ** (1.0_dp / shape)
         return
      end if

      a = shape
      d = a - 1.0_dp / 3.0_dp
      c = 1.0_dp / sqrt(9.0_dp * d)
      do
         z = random_standard_normal()
         v = 1.0_dp + c * z
         if (v <= 0.0_dp) cycle
         v = v * v * v
         call random_number(u)
         if (u < 1.0_dp - 0.0331_dp * z ** 4) exit
         if (log(max(u, tiny(1.0_dp))) < 0.5_dp * z * z + d * (1.0_dp - v + log(v))) exit
      end do
      x = scl * d * v
   end function random_gamma

   function random_poisson(lambda) result(k)
      ! Exact inversion for small means and PTRS transformed rejection for
      ! larger means (Hormann, 1993).
      real(dp), intent(in) :: lambda
      integer :: k
      real(dp) :: a, b, inv_alpha, log_lambda, p, s, u, us, v, vr

      if (lambda <= 0.0_dp) then
         k = 0
         return
      end if

      if (lambda < 30.0_dp) then
         p = exp(-lambda)
         s = p
         call random_number(u)
         k = 0
         do while (u > s)
            k = k + 1
            p = p * lambda / real(k, dp)
            s = s + p
         end do
         return
      end if

      s = sqrt(lambda)
      b = 0.931_dp + 2.53_dp * s
      a = -0.059_dp + 0.02483_dp * b
      inv_alpha = 1.1239_dp + 1.1328_dp / (b - 3.4_dp)
      vr = 0.9277_dp - 3.6224_dp / (b - 2.0_dp)
      log_lambda = log(lambda)

      do
         call random_number(u)
         call random_number(v)
         u = u - 0.5_dp
         us = 0.5_dp - abs(u)
         if (us <= tiny(1.0_dp)) cycle
         k = floor((2.0_dp * a / us + b) * u + lambda + 0.43_dp)
         if (us >= 0.07_dp .and. v <= vr) return
         if (k < 0) cycle
         if (us < 0.013_dp .and. v > us) cycle
         if (log(v * inv_alpha / (a / (us * us) + b)) <= &
             -lambda + real(k, dp) * log_lambda - log_gamma(real(k + 1, dp))) return
      end do
   end function random_poisson

   subroutine correlation_matrix(x, corr, status, message)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: corr(:,:)
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      integer :: i, j, n, p
      real(dp), allocatable :: mean_x(:), sd_x(:)
      real(dp) :: denom

      n = size(x, 1)
      p = size(x, 2)
      status = 0
      message = ''
      if (n < 2 .or. p < 1) then
         status = 1
         message = 'correlation_matrix requires at least two rows and one column'
         allocate(corr(0, 0))
         return
      end if

      allocate(corr(p, p), mean_x(p), sd_x(p))
      mean_x = sum(x, dim=1) / real(n, dp)
      do j = 1, p
         sd_x(j) = sqrt(sum((x(:, j) - mean_x(j)) ** 2) / real(n - 1, dp))
      end do

      corr = 0.0_dp
      do i = 1, p
         corr(i, i) = 1.0_dp
         do j = i + 1, p
            denom = real(n - 1, dp) * sd_x(i) * sd_x(j)
            if (denom <= tiny(1.0_dp)) then
               status = 2
               message = 'a scenario-factor column has zero variance'
               corr(i, j) = 0.0_dp
            else
               corr(i, j) = sum((x(:, i) - mean_x(i)) * &
                                (x(:, j) - mean_x(j))) / denom
            end if
            corr(j, i) = corr(i, j)
         end do
      end do
   end subroutine correlation_matrix

   pure function quadratic_form(x, a) result(value)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: value

      value = dot_product(x, matmul(a, x))
   end function quadratic_form

end module gcpm_math
