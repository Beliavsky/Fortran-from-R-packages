program test_negbin_change
   use mcmcpack
   implicit none
   integer,parameter::n=12,ns=2,k=2
   integer::y(n),i
   real(dp)::x(n,k),beta(ns,k),rho(ns),p(ns,ns),a0(ns,ns),b0(k),bp(k,k),step(ns)
   type(change_result)::r
   x(:,1)=1.0_dp;do i=1,n;x(i,2)=real(i-6,dp)/6.0_dp;end do
   y=[0,1,0,1,1,2,2,3,4,2,5,4]
   beta=0.0_dp;rho=[2.0_dp,3.0_dp];p=0.0_dp;p(1,:)=[0.85_dp,0.15_dp];p(2,2)=1.0_dp
   a0=0.0_dp;a0(1,:)=[8.0_dp,1.0_dp];a0(2,2)=1.0_dp
   b0=0.0_dp;bp=0.0_dp;bp(1,1)=0.1_dp;bp(2,2)=0.1_dp;step=0.5_dp
   call set_seed(1423)
   r=mcmc_negbin_change(y,x,beta,rho,p,b0,bp,a0,2.0_dp,2.0_dp,1.0_dp,step,20,40,2,0.25_dp)
   if(r%status/=0)error stop 1
   if(any(shape(r%draws)/=[20,ns*k+ns+ns*ns]))error stop 2
   if(any(r%draws(:,ns*k+1:ns*k+ns)<=0.0_dp))error stop 3
   if(any(r%states<1).or.any(r%states>ns))error stop 4
   print '(a)','test_negbin_change: PASS'
end program test_negbin_change
