program test_negbin
   use mcmcpack
   implicit none
   integer,parameter::n=12,k=2
   integer::y(n),i
   real(dp)::x(n,k),b(k),b0(k),bp(k,k)
   type(mcmc_result)::r
   call set_seed(13579)
   y=[0,1,0,2,1,3,0,1,4,2,1,0]
   do i=1,n;x(i,1)=1.0_dp;x(i,2)=real(i-6,dp)/5.0_dp;end do
   b=[0.0_dp,0.0_dp];b0=0.0_dp;bp=0.0_dp;bp(1,1)=0.1_dp;bp(2,2)=0.1_dp
   r=mcmc_negbin(y,x,b,2.0_dp,b0,bp,2.0_dp,2.0_dp,1.0_dp,1.0_dp,0.2_dp,20,40,2)
   if(r%status/=0.or.size(r%draws,1)/=20.or.size(r%draws,2)/=3) error stop 'negbin'
   if(any(r%draws(:,3)<=0.0_dp)) error stop 'negbin rho'
   print '(a)','test_negbin: PASS'
end program
