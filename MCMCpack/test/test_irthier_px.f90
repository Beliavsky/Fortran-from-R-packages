program test_irthier_px
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : set_seed
   use mcmcpack_irthier, only : irthier_result,mcmc_irt_hier1d
   implicit none
   integer::x(4,3)
   real(dp)::th(4),eta(3,2),ab0(2),abp(2,2),xj(4,2),bst(2),b0(2),bp(2,2)
   type(irthier_result)::r
   call set_seed(314159)
   x=reshape([1,0,1,0, 0,1,1,0, 1,1,0,1],[4,3])
   th=[-0.5_dp,0.2_dp,0.6_dp,-0.1_dp];eta=0.0_dp;eta(:,2)=1.0_dp;ab0=0.0_dp;abp=0.0_dp;abp(1,1)=1.0_dp;abp(2,2)=1.0_dp
   xj(:,1)=1.0_dp;xj(:,2)=[-1.0_dp,-0.3_dp,0.4_dp,1.0_dp];bst=0.0_dp;b0=0.0_dp;bp=0.0_dp;bp(1,1)=1.0_dp;bp(2,2)=1.0_dp
   r=mcmc_irt_hier1d(x,th,eta,ab0,abp,xj,bst,b0,bp,2.0_dp,2.0_dp,5,10,1,.true.,.true.,.true.,2.0_dp,2.0_dp,.true.)
   if(r%status/=0.or.size(r%draws,1)/=10.or..not.(r%log_marginal==r%log_marginal))error stop 1
   print *,'test_irthier_px: PASS'
end program
