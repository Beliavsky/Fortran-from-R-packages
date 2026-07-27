! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_random
   use sde_kinds, only : dp, i64, pi, tiny_dp
   implicit none
   private

   public :: seed_rng
   public :: random_uniform
   public :: random_normal
   public :: random_gamma
   public :: random_chi_square
   public :: random_poisson
   public :: random_noncentral_chi_square

   logical, save :: has_spare_normal = .false.
   real(dp), save :: spare_normal = 0.0_dp

contains

   subroutine seed_rng(seed)
      integer(i64), intent(in) :: seed
      integer, allocatable :: values(:)
      integer :: n, i
      integer(i64) :: work

      call random_seed(size=n)
      allocate(values(n))
      do i = 1, n
         work = modulo(abs(seed)+104729_i64*int(i, i64)+8191_i64*int(i*i, i64), &
            int(huge(1), i64)-1_i64)
         values(i) = int(work+1_i64)
      end do
      call random_seed(put=values)
      has_spare_normal = .false.
   end subroutine seed_rng

   function random_uniform() result(value)
      real(dp) :: value

      call random_number(value)
      value = max(tiny_dp, min(1.0_dp-tiny_dp, value))
   end function random_uniform

   function random_normal(mean, sd) result(value)
      real(dp), intent(in), optional :: mean
      real(dp), intent(in), optional :: sd
      real(dp) :: value
      real(dp) :: location, scale, u1, u2, radius

      location = 0.0_dp
      scale = 1.0_dp
      if (present(mean)) location = mean
      if (present(sd)) scale = sd
      if (scale < 0.0_dp) error stop "random_normal: sd must be nonnegative"

      if (has_spare_normal) then
         value = spare_normal
         has_spare_normal = .false.
      else
         u1 = random_uniform()
         u2 = random_uniform()
         radius = sqrt(-2.0_dp*log(u1))
         value = radius*cos(2.0_dp*pi*u2)
         spare_normal = radius*sin(2.0_dp*pi*u2)
         has_spare_normal = .true.
      end if
      value = location+scale*value
   end function random_normal

   recursive function random_gamma(shape, scale) result(value)
      real(dp), intent(in) :: shape
      real(dp), intent(in), optional :: scale
      real(dp) :: value
      real(dp) :: scale_value, d, c, x, v, u

      scale_value = 1.0_dp
      if (present(scale)) scale_value = scale
      if (shape <= 0.0_dp .or. scale_value <= 0.0_dp) then
         error stop "random_gamma: shape and scale must be positive"
      end if

      if (shape < 1.0_dp) then
         value = random_gamma(shape+1.0_dp)*random_uniform()**(1.0_dp/shape)
         value = scale_value*value
         return
      end if

      d = shape-1.0_dp/3.0_dp
      c = 1.0_dp/sqrt(9.0_dp*d)
      do
         do
            x = random_normal()
            v = 1.0_dp+c*x
            if (v > 0.0_dp) exit
         end do
         v = v*v*v
         u = random_uniform()
         if (u < 1.0_dp-0.0331_dp*x**4) exit
         if (log(u) < 0.5_dp*x*x+d*(1.0_dp-v+log(v))) exit
      end do
      value = scale_value*d*v
   end function random_gamma

   function random_chi_square(df) result(value)
      real(dp), intent(in) :: df
      real(dp) :: value

      if (df <= 0.0_dp) error stop "random_chi_square: df must be positive"
      value = random_gamma(0.5_dp*df, 2.0_dp)
   end function random_chi_square

   function random_poisson(lambda) result(value)
      real(dp), intent(in) :: lambda
      integer :: value
      real(dp) :: product, u, v, us, b, a, inv_alpha, vr, log_lambda
      integer :: candidate

      if (lambda < 0.0_dp) error stop "random_poisson: lambda must be nonnegative"
      if (lambda <= 0.0_dp) then
         value = 0
         return
      end if

      if (lambda < 30.0_dp) then
         product = 1.0_dp
         value = -1
         do
            value = value+1
            product = product*random_uniform()
            if (product <= exp(-lambda)) exit
         end do
         return
      end if

      ! Hormann's PTRS transformed rejection algorithm.
      b = 0.931_dp+2.53_dp*sqrt(lambda)
      a = -0.059_dp+0.02483_dp*b
      inv_alpha = 1.1239_dp+1.1328_dp/(b-3.4_dp)
      vr = 0.9277_dp-3.6224_dp/(b-2.0_dp)
      log_lambda = log(lambda)
      do
         u = random_uniform()-0.5_dp
         v = random_uniform()
         us = 0.5_dp-abs(u)
         candidate = floor((2.0_dp*a/us+b)*u+lambda+0.43_dp)
         if (us >= 0.07_dp .and. v <= vr .and. candidate >= 0) then
            value = candidate
            return
         end if
         if (candidate < 0 .or. (us < 0.013_dp .and. v > us)) cycle
         if (log(v*inv_alpha/(a/(us*us)+b)) <= &
             -lambda+real(candidate, dp)*log_lambda-log_gamma(real(candidate+1, dp))) then
            value = candidate
            return
         end if
      end do
   end function random_poisson

   function random_noncentral_chi_square(df, ncp) result(value)
      real(dp), intent(in) :: df
      real(dp), intent(in) :: ncp
      real(dp) :: value
      integer :: k

      if (df <= 0.0_dp .or. ncp < 0.0_dp) then
         error stop "random_noncentral_chi_square: invalid df or ncp"
      end if
      k = random_poisson(0.5_dp*ncp)
      value = random_chi_square(df+2.0_dp*real(k, dp))
   end function random_noncentral_chi_square

end module sde_random
