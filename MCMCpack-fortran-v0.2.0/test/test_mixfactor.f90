program test_mixfactor
   use mcmcpack
   implicit none
   integer,parameter::n=12,k=3,d=2,gdim=4
   real(dp)::x(n,k),lam(k,d),gam(gdim,k),psi(k),eq(k,d),ineq(k,d),pm(k,d),pp(k,d),a0(k),b0(k),tune(k)
   integer::ncat(k),i
   type(mixfactor_result)::r
   do i=1,n;x(i,1)=0.2_dp*real(i-6,dp);end do
   x(:,2)=real([1,2,3,1,2,3,1,2,3,1,2,3],dp);x(:,3)=real([1,2,1,2,1,2,1,2,1,2,1,2],dp)
   ncat=[0,3,2];lam=0.0_dp;lam(:,2)=[0.8_dp,0.7_dp,-0.6_dp];psi=[1.0_dp,1.0_dp,1.0_dp]
   gam=300.0_dp;do i=1,k;gam(1,i)=-300.0_dp;gam(2,i)=0.0_dp;end do;gam(3,2)=0.8_dp;gam(4,2)=300.0_dp
   eq=-999.0_dp;ineq=0.0_dp;eq(1,2)=0.8_dp;pm=0.0_dp;pp=1.0_dp;a0=2.0_dp;b0=2.0_dp;tune=[0.0_dp,0.08_dp,0.0_dp]
   call set_seed(1191)
   r=mcmc_mixfactanal(x,ncat,lam,gam,psi,eq,ineq,pm,pp,a0,b0,tune,20,40,2,.true.,.true.)
   if(r%status/=0)error stop 1
   if(any(shape(r%draws)/=[20,k*d+gdim*k+n*d+k]))error stop 2
   if(any(r%draws(:,size(r%draws,2)-k+1:size(r%draws,2))<=0.0_dp))error stop 3
   print '(a)','test_mixfactor: PASS'
end program
