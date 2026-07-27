! Part of the experimental modern Fortran translation of tseries 0.10-62.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original tseries authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only

module tseries_random
   use tseries_kinds, only : dp
   implicit none
   private

   public :: seed_random
   public :: random_normal
   public :: random_permutation
   public :: random_integer
   public :: random_exponential

contains

   subroutine seed_random(seed)
      integer, intent(in) :: seed
      integer, allocatable :: put(:)
      integer :: n,i
      call random_seed(size=n)
      allocate(put(n))
      do i=1,n
         put(i)=modulo(seed+104729*i,huge(1)-1)+1
      end do
      call random_seed(put=put)
   end subroutine seed_random

   real(dp) function random_normal() result(z)
      real(dp) :: u1,u2
      real(dp), parameter :: two_pi=2.0_dp*acos(-1.0_dp)
      call random_number(u1)
      call random_number(u2)
      u1=max(u1,tiny(1.0_dp))
      z=sqrt(-2.0_dp*log(u1))*cos(two_pi*u2)
   end function random_normal

   real(dp) function random_exponential() result(x)
      real(dp) :: u
      call random_number(u)
      x=-log(max(u,tiny(1.0_dp)))
   end function random_exponential

   integer function random_integer(low,high) result(value)
      integer, intent(in) :: low,high
      real(dp) :: u
      if(high<low) then
         value=low; return
      end if
      call random_number(u)
      value=low+int(u*real(high-low+1,dp))
      value=min(value,high)
   end function random_integer

   subroutine random_permutation(n,p)
      integer, intent(in) :: n
      integer, intent(out) :: p(:)
      integer :: i,j,tmp
      if(size(p)/=n) return
      p=[(i,i=1,n)]
      do i=n,2,-1
         j=random_integer(1,i)
         tmp=p(i); p(i)=p(j); p(j)=tmp
      end do
   end subroutine random_permutation

end module tseries_random
