! Part of the experimental modern Fortran translation of fGarch 4052.93.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original fGarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

module fgarch_rng
   use fgarch_kinds, only : dp
   use fgarch_math, only : pi
   implicit none
   private

   public :: seed_rng, random_normal, random_gamma, random_student_t

contains

   subroutine seed_rng(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)

      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(seed + 104729*i + 37*i*i, huge(1)-1)
         if (put(i) <= 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine seed_rng

   function random_normal() result(z)
      real(dp) :: z
      real(dp) :: u1, u2

      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
   end function random_normal

   recursive function random_gamma(shape) result(x)
      real(dp), intent(in) :: shape
      real(dp) :: x
      real(dp) :: d, c, z, u, v

      if (shape <= 0.0_dp) then
         x = 0.0_dp
         return
      end if
      if (shape < 1.0_dp) then
         call random_number(u)
         x = random_gamma(shape+1.0_dp)*u**(1.0_dp/shape)
         return
      end if

      d = shape - 1.0_dp/3.0_dp
      c = 1.0_dp/sqrt(9.0_dp*d)
      do
         do
            z = random_normal()
            v = 1.0_dp + c*z
            if (v > 0.0_dp) exit
         end do
         v = v*v*v
         call random_number(u)
         if (u < 1.0_dp-0.0331_dp*z**4) exit
         if (log(max(u,tiny(1.0_dp))) < 0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
      end do
      x = d*v
   end function random_gamma

   function random_student_t(nu) result(x)
      real(dp), intent(in) :: nu
      real(dp) :: x, chi2

      chi2 = 2.0_dp*random_gamma(0.5_dp*nu)
      x = random_normal()/sqrt(chi2/nu)
   end function random_student_t

end module fgarch_rng
