! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! Derived from parma 1.7, Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
module parma_rng
   use parma_kinds, only: dp, pi
   implicit none
   private
   public :: seed_rng, random_normal, random_normals, random_rademacher

contains

   subroutine seed_rng(seed)
      integer, intent(in) :: seed
      integer :: i, n
      integer, allocatable :: seeds(:)

      call random_seed(size=n)
      allocate(seeds(n))
      do i = 1, n
         seeds(i) = modulo(seed + 104729 * i, huge(1) - 1)
         if (seeds(i) <= 0) seeds(i) = i
      end do
      call random_seed(put=seeds)
   end subroutine seed_rng

   function random_normal() result(z)
      real(dp) :: z
      real(dp) :: u1, u2

      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
   end function random_normal

   subroutine random_normals(z)
      real(dp), intent(out) :: z(:)
      integer :: i
      do i = 1, size(z)
         z(i) = random_normal()
      end do
   end subroutine random_normals

   function random_rademacher() result(x)
      real(dp) :: x
      real(dp) :: u
      call random_number(u)
      if (u < 0.5_dp) then
         x = -1.0_dp
      else
         x = 1.0_dp
      end if
   end function random_rademacher

end module parma_rng
