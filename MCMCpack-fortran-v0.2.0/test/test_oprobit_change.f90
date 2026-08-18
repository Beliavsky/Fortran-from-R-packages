program test_oprobit_change
   use mcmcpack
   implicit none
   integer,parameter::n=12,ns=2,k=2,gk=4
   integer::y(n),i
   real(dp)::x(n,k),b(ns,k),bl(ns,k),gam(ns,gk),p(ns,ns),b0(k),bp(k,k),a0(ns,ns),tune(ns)
   type(oprobit_change_result)::r
   x(:,1)=1.0_dp;do i=1,n;x(i,2)=real(i-6,dp)/5.0_dp;end do
   y=[1,1,2,1,2,2,2,3,2,3,3,3];b=0.0_dp;bl=0.0_dp
   gam=0.0_dp;gam(:,1)=-300.0_dp;gam(:,2)=0.0_dp;gam(:,3)=0.8_dp;gam(:,4)=300.0_dp
   p=0.0_dp;p(1,:)=[0.85_dp,0.15_dp];p(2,2)=1.0_dp
   b0=0.0_dp;bp=0.0_dp;bp(1,1)=0.1_dp;bp(2,2)=0.1_dp
   a0=0.0_dp;a0(1,:)=[8.0_dp,1.0_dp];a0(2,2)=1.0_dp;tune=0.08_dp
   call set_seed(7511)
   r=mcmc_oprobit_change(y,x,b,bl,gam,p,1.0_dp,b0,bp,a0,tune,20,40,2,.false.)
   if(r%status/=0)error stop 1
   if(any(shape(r%draws)/=[20,2*ns*k+ns*gk+ns*ns]))error stop 2
   if(any(r%states<1).or.any(r%states>ns))error stop 3
   if(any(r%gamma_accept_rate<0.0_dp).or.any(r%gamma_accept_rate>1.0_dp))error stop 4
   print '(a)','test_oprobit_change: PASS'
end program
