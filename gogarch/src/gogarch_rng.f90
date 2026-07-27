! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module gogarch_rng
   use gogarch_kinds, only : dp
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)
   logical :: has_spare = .false.
   real(dp) :: spare = 0.0_dp
   public :: seed_rng, random_uniform, random_normal, random_gamma
   public :: random_student_t, fill_random_normal
contains
   subroutine seed_rng(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)
      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(seed+104729*i,huge(1)-1)+1
      end do
      call random_seed(put=put)
      has_spare = .false.
   end subroutine seed_rng

   function random_uniform() result(u)
      real(dp) :: u
      call random_number(u)
      u = max(tiny(1.0_dp),min(1.0_dp-tiny(1.0_dp),u))
   end function random_uniform

   function random_normal() result(z)
      real(dp) :: z, u1, u2, radius
      if (has_spare) then
         z = spare
         has_spare = .false.
         return
      end if
      u1 = random_uniform()
      u2 = random_uniform()
      radius = sqrt(-2.0_dp*log(u1))
      z = radius*cos(2.0_dp*pi*u2)
      spare = radius*sin(2.0_dp*pi*u2)
      has_spare = .true.
   end function random_normal

   recursive function random_gamma(shape) result(x)
      real(dp), intent(in) :: shape
      real(dp) :: x, d, c, z, v, u
      if (shape <= 0.0_dp) error stop 'random_gamma: shape must be positive'
      if (shape < 1.0_dp) then
         x = random_gamma(shape+1.0_dp)*random_uniform()**(1.0_dp/shape)
         return
      end if
      d = shape-1.0_dp/3.0_dp
      c = 1.0_dp/sqrt(9.0_dp*d)
      do
         do
            z = random_normal()
            v = 1.0_dp+c*z
            if (v > 0.0_dp) exit
         end do
         v = v**3
         u = random_uniform()
         if (u < 1.0_dp-0.0331_dp*z**4) exit
         if (log(u) < 0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
      end do
      x = d*v
   end function random_gamma

   function random_student_t(df) result(x)
      real(dp), intent(in) :: df
      real(dp) :: x
      if (df <= 0.0_dp) error stop 'random_student_t: df must be positive'
      x = random_normal()/sqrt(2.0_dp*random_gamma(0.5_dp*df)/df)
   end function random_student_t

   subroutine fill_random_normal(x)
      real(dp), intent(out) :: x(..)
      integer :: i, j
      select rank (x)
      rank (1)
         do i = 1, size(x)
            x(i) = random_normal()
         end do
      rank (2)
         do j = 1, size(x,2)
            do i = 1, size(x,1)
               x(i,j) = random_normal()
            end do
         end do
      end select
   end subroutine fill_random_normal
end module gogarch_rng
