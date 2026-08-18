program test_irthier
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : set_seed
   use mcmcpack_irthier, only : irthier_result,mcmc_irt_hier1d
   implicit none
   integer :: x(6,3)
   real(dp) :: th(6),eta(3,2),ab0(2),abp(2,2),xj(6,2),be(2),b0(2),bp(2,2)
   type(irthier_result) :: r
   integer :: i
   x=reshape([1,1,0,0,1,0, 0,1,1,0,1,0, 1,0,1,1,0,0],[6,3])
   th=[(-0.75_dp+0.3_dp*real(i-1,dp),i=1,6)]
   eta(:,1)=0.0_dp;eta(:,2)=1.0_dp;ab0=[0.0_dp,1.0_dp];abp=0.0_dp;abp(1,1)=0.1_dp;abp(2,2)=0.1_dp
   xj(:,1)=1.0_dp;xj(:,2)=[(-1.0_dp+0.4_dp*real(i-1,dp),i=1,6)]
   be=0.0_dp;b0=0.0_dp;bp=0.0_dp;bp(1,1)=0.1_dp;bp(2,2)=0.1_dp
   call set_seed(8181)
   r=mcmc_irt_hier1d(x,th,eta,ab0,abp,xj,be,b0,bp,2.0_dp,2.0_dp,10,20,2)
   if(r%status/=0) error stop 'irthier status'
   if(any(shape(r%draws)/=[10,15])) error stop 'irthier shape'
   if(any(r%draws(:,15)<=0.0_dp)) error stop 'irthier sigma2'
   print *,'test_irthier: PASS'
end program
