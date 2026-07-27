! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2008-2021 David Ardia
! Modern Fortran translation of computational routines from bayesGARCH.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
module bayesgarch_math
   use bayesgarch_kinds, only : dp
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)

   public :: digamma
   public :: inverse_2x2
   public :: cholesky_2x2
   public :: log_normal_density
   public :: log_mvn2_density
   public :: draw_mvn2
   public :: student_log_density

contains

   pure function digamma(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      real(dp) :: xx
      real(dp) :: inv
      real(dp) :: inv2

      if (x <= 0.0_dp) then
         value = huge(1.0_dp)
         return
      end if

      value = 0.0_dp
      xx = x
      do while (xx < 8.0_dp)
         value = value - 1.0_dp / xx
         xx = xx + 1.0_dp
      end do

      inv = 1.0_dp / xx
      inv2 = inv * inv
      value = value + log(xx) - 0.5_dp * inv - inv2 * (1.0_dp / 12.0_dp - &
         inv2 * (1.0_dp / 120.0_dp - inv2 * (1.0_dp / 252.0_dp - &
         inv2 * (1.0_dp / 240.0_dp - inv2 * (5.0_dp / 660.0_dp)))))
   end function digamma

   pure subroutine inverse_2x2(a, ainv, determinant, ok)
      real(dp), intent(in) :: a(2, 2)
      real(dp), intent(out) :: ainv(2, 2)
      real(dp), intent(out) :: determinant
      logical, intent(out) :: ok

      determinant = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
      ok = determinant > tiny(1.0_dp)
      if (.not. ok) then
         ainv = 0.0_dp
         return
      end if
      ainv(1, 1) = a(2, 2) / determinant
      ainv(1, 2) = -a(1, 2) / determinant
      ainv(2, 1) = -a(2, 1) / determinant
      ainv(2, 2) = a(1, 1) / determinant
   end subroutine inverse_2x2

   pure subroutine cholesky_2x2(a, l, ok)
      real(dp), intent(in) :: a(2, 2)
      real(dp), intent(out) :: l(2, 2)
      logical, intent(out) :: ok
      real(dp) :: remainder

      l = 0.0_dp
      ok = a(1, 1) > 0.0_dp
      if (.not. ok) return
      l(1, 1) = sqrt(a(1, 1))
      l(2, 1) = a(2, 1) / l(1, 1)
      remainder = a(2, 2) - l(2, 1)**2
      ok = remainder > 0.0_dp
      if (.not. ok) then
         l = 0.0_dp
         return
      end if
      l(2, 2) = sqrt(remainder)
   end subroutine cholesky_2x2

   pure function log_normal_density(x, mean, variance) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: mean
      real(dp), intent(in) :: variance
      real(dp) :: value

      if (variance <= 0.0_dp) then
         value = -huge(1.0_dp)
      else
         value = -0.5_dp * (log(2.0_dp * pi * variance) + (x - mean)**2 / variance)
      end if
   end function log_normal_density

   pure function log_mvn2_density(x, mean, covariance) result(value)
      real(dp), intent(in) :: x(2)
      real(dp), intent(in) :: mean(2)
      real(dp), intent(in) :: covariance(2, 2)
      real(dp) :: value
      real(dp) :: inv(2, 2)
      real(dp) :: determinant
      real(dp) :: d(2)
      logical :: ok

      call inverse_2x2(covariance, inv, determinant, ok)
      if (.not. ok) then
         value = -huge(1.0_dp)
         return
      end if
      d = x - mean
      value = -log(2.0_dp * pi) - 0.5_dp * log(determinant) - &
         0.5_dp * dot_product(d, matmul(inv, d))
   end function log_mvn2_density

   function draw_mvn2(mean, covariance) result(x)
      use bayesgarch_rng, only : random_normal
      real(dp), intent(in) :: mean(2)
      real(dp), intent(in) :: covariance(2, 2)
      real(dp) :: x(2)
      real(dp) :: l(2, 2)
      real(dp) :: z(2)
      logical :: ok

      call cholesky_2x2(covariance, l, ok)
      if (.not. ok) error stop "draw_mvn2: covariance is not positive definite"
      z = [random_normal(), random_normal()]
      x = mean + matmul(l, z)
   end function draw_mvn2

   pure function student_log_density(x, variance, nu) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: variance
      real(dp), intent(in) :: nu
      real(dp) :: value

      if (variance <= 0.0_dp .or. nu <= 2.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      value = log_gamma(0.5_dp * (nu + 1.0_dp)) - log_gamma(0.5_dp * nu) - &
         0.5_dp * log(pi * (nu - 2.0_dp) * variance) - &
         0.5_dp * (nu + 1.0_dp) * log(1.0_dp + x * x / ((nu - 2.0_dp) * variance))
   end function student_log_density

end module bayesgarch_math
