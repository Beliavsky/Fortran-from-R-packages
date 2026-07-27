! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module msgarch_mapping
   use msgarch_kinds, only : dp
   implicit none
   private
   public :: bounded_map, bounded_unmap, simplex_mapping, simplex_unmapping
contains
   pure function bounded_map(x,lower,upper) result(value)
      real(dp),intent(in)::x,lower,upper
      real(dp)::value
      if(x>=0.0_dp)then
         value=lower+(upper-lower)/(1.0_dp+exp(-x))
      else
         value=lower+(upper-lower)*exp(x)/(1.0_dp+exp(x))
      end if
   end function bounded_map

   pure function bounded_unmap(value,lower,upper) result(x)
      real(dp),intent(in)::value,lower,upper
      real(dp)::x,p
      p=max(1.0e-14_dp,min(1.0_dp-1.0e-14_dp,(value-lower)/(upper-lower)))
      x=log(p)-log(1.0_dp-p)
   end function bounded_unmap

   function simplex_mapping(phi) result(probability)
      real(dp),intent(in)::phi(:)
      real(dp),allocatable::probability(:)
      real(dp)::remaining,p
      integer::i,k
      k=size(phi)+1;allocate(probability(k));remaining=1.0_dp
      do i=1,k-1
         if(phi(i)>=0.0_dp)then;p=1.0_dp/(1.0_dp+exp(-phi(i)));else;p=exp(phi(i))/(1.0_dp+exp(phi(i)));end if
         probability(i)=remaining*p;remaining=remaining*(1.0_dp-p)
      end do
      probability(k)=remaining
   end function simplex_mapping

   function simplex_unmapping(probability) result(phi)
      real(dp),intent(in)::probability(:)
      real(dp),allocatable::phi(:)
      real(dp)::remaining,p
      integer::i,k
      k=size(probability);allocate(phi(k-1));remaining=1.0_dp
      do i=1,k-1
         p=max(1.0e-14_dp,min(1.0_dp-1.0e-14_dp,probability(i)/remaining))
         phi(i)=log(p)-log(1.0_dp-p);remaining=remaining-probability(i)
      end do
   end function simplex_unmapping
end module msgarch_mapping
