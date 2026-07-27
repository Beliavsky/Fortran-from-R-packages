! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2008-2021 David Ardia
! Modern Fortran translation of computational routines from bayesGARCH.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
module bayesgarch_rng
   use bayesgarch_kinds, only : dp
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)

   public :: seed_rng
   public :: random_uniform
   public :: random_normal
   public :: random_exponential
   public :: random_gamma
   public :: random_standardized_student

contains

   subroutine seed_rng(seed)
      integer, intent(in) :: seed
      integer :: i
      integer :: n
      integer, allocatable :: put(:)
      integer(kind=8) :: value

      call random_seed(size=n)
      allocate(put(n))
      value = int(seed, kind=8)
      do i = 1, n
         value = modulo(1103515245_8 * value + 12345_8 + 104729_8 * int(i, 8), 2147483647_8)
         if (value == 0_8) value = int(i, 8)
         put(i) = int(value)
      end do
      call random_seed(put=put)
   end subroutine seed_rng

   function random_uniform() result(x)
      real(dp) :: x

      call random_number(x)
      if (x <= 0.0_dp) x = tiny(1.0_dp)
      if (x >= 1.0_dp) x = 1.0_dp - epsilon(1.0_dp)
   end function random_uniform

   function random_normal() result(x)
      real(dp) :: x
      real(dp) :: u1
      real(dp) :: u2

      u1 = random_uniform()
      u2 = random_uniform()
      x = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
   end function random_normal

   function random_exponential(rate) result(x)
      real(dp), intent(in) :: rate
      real(dp) :: x

      if (rate <= 0.0_dp) error stop "random_exponential: rate must be positive"
      x = -log(random_uniform()) / rate
   end function random_exponential

   recursive function random_gamma(shape, rate) result(x)
      real(dp), intent(in) :: shape
      real(dp), intent(in) :: rate
      real(dp) :: x
      real(dp) :: c
      real(dp) :: d
      real(dp) :: u
      real(dp) :: v
      real(dp) :: z

      if (shape <= 0.0_dp) error stop "random_gamma: shape must be positive"
      if (rate <= 0.0_dp) error stop "random_gamma: rate must be positive"

      if (shape < 1.0_dp) then
         x = random_gamma(shape + 1.0_dp, rate) * random_uniform()**(1.0_dp / shape)
         return
      end if

      d = shape - 1.0_dp / 3.0_dp
      c = 1.0_dp / sqrt(9.0_dp * d)
      do
         do
            z = random_normal()
            v = 1.0_dp + c * z
            if (v > 0.0_dp) exit
         end do
         v = v**3
         u = random_uniform()
         if (u < 1.0_dp - 0.0331_dp * z**4) exit
         if (log(u) < 0.5_dp * z * z + d * (1.0_dp - v + log(v))) exit
      end do
      x = d * v / rate
   end function random_gamma

   function random_standardized_student(nu) result(x)
      real(dp), intent(in) :: nu
      real(dp) :: x
      real(dp) :: chi2

      if (nu <= 2.0_dp) error stop "random_standardized_student: nu must exceed 2"
      chi2 = 2.0_dp * random_gamma(0.5_dp * nu, 1.0_dp)
      x = random_normal() * sqrt((nu - 2.0_dp) / chi2)
   end function random_standardized_student

end module bayesgarch_rng
