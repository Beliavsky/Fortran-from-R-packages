! SPDX-License-Identifier: GPL-3.0-only
! Numerical translation of fHMM 1.4.3.
module fhmm_calendar
   use fhmm_kinds, only: dp
   use fhmm_math, only: seed_rng
   implicit none
   private
   public :: compute_chunk_lengths
contains
   function compute_chunk_lengths(coarse_horizon,fine_horizon,period,seed) result(lengths)
      integer,intent(in)::coarse_horizon,fine_horizon
      character(len=*),intent(in)::period
      integer,intent(in),optional::seed
      integer,allocatable::lengths(:)
      integer::size_hint,t,k
      real(dp)::u,c,total
      real(dp),allocatable::weights(:)
      allocate(lengths(coarse_horizon))
      if(fine_horizon>0)then;lengths=fine_horizon;return;end if
      select case(trim(period))
      case('w')
         size_hint=5
      case('m')
         size_hint=25
      case('q')
         size_hint=70
      case('y')
         size_hint=260
      case default
         size_hint=25
      end select
      if(present(seed))call seed_rng(seed)
      allocate(weights(size_hint));total=0.0_dp
      do k=1,size_hint
         weights(k)=binomial_probability(k,size_hint,0.9_dp);total=total+weights(k)
      end do
      weights=weights/total
      do t=1,coarse_horizon
         call random_number(u);c=0.0_dp;lengths(t)=size_hint
         do k=1,size_hint;c=c+weights(k);if(u<=c)then;lengths(t)=k;exit;end if;end do
      end do
   end function compute_chunk_lengths
   pure real(dp) function binomial_probability(k,n,p) result(v)
      integer,intent(in)::k,n
      real(dp),intent(in)::p
      v=exp(log_gamma(real(n+1,dp))-log_gamma(real(k+1,dp))-log_gamma(real(n-k+1,dp))+ &
         real(k,dp)*log(p)+real(n-k,dp)*log(1.0_dp-p))
   end function binomial_probability
end module fhmm_calendar
