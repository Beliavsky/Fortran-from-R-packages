! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from OptionPricing 0.1.2 by Wolfgang Hormann and Kemal Dingec.
module optionpricing_random
   use optionpricing_kinds, only : dp, pi
   implicit none
   private
   public :: seed_rng, random_normal, fill_normal
contains
   subroutine seed_rng(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)
      call random_seed(size=n)
      allocate(put(n))
      do i=1,n
         put(i) = modulo(seed + 104729*i + 37*i*i, huge(1)-1)
         if (put(i) <= 0) put(i) = i+seed+1
      end do
      call random_seed(put=put)
   end subroutine seed_rng

   function random_normal() result(z)
      real(dp) :: z, u1, u2
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
   end function random_normal

   subroutine fill_normal(x)
      real(dp), intent(out) :: x(..)
      integer :: i, j
      select rank(x)
      rank(1)
         do i=1,size(x)
            x(i) = random_normal()
         end do
      rank(2)
         do j=1,size(x,2)
            do i=1,size(x,1)
               x(i,j) = random_normal()
            end do
         end do
      end select
   end subroutine fill_normal
end module optionpricing_random
