program test_ordfactor
   use mcmcpack
   implicit none
   integer,parameter::n=12,k=3,d=2,gdim=4
   integer::x(n,k),ncat(k),i
   real(dp)::lam(k,d),gam(gdim,k),eq(k,d),ineq(k,d),pm(k,d),pp(k,d),tune(k)
   type(ordfactor_result)::a
   x(:,1)=[1,1,2,2,1,2,1,2,2,1,2,1]
   x(:,2)=[1,2,3,1,2,3,1,2,3,1,2,3]
   x(:,3)=[2,1,2,1,2,1,2,1,2,1,2,1]
   ncat=[2,3,2];lam=0.0_dp;lam(:,2)=[1.0_dp,0.7_dp,-0.5_dp]
   gam=300.0_dp
   do i=1,k;gam(1,i)=-300.0_dp;gam(2,i)=0.0_dp;end do
   gam(3,1)=300.0_dp;gam(3,2)=0.8_dp;gam(4,2)=300.0_dp;gam(3,3)=300.0_dp
   eq=-999.0_dp;ineq=0.0_dp;eq(1,2)=1.0_dp;ineq(2,2)=1.0_dp
   pm=0.0_dp;pp=1.0_dp;tune=[0.0_dp,0.08_dp,0.0_dp]
   call set_seed(9124)
   a=mcmc_ordfactanal(x,lam,gam,ncat,eq,ineq,pm,pp,tune,20,40,2,.true.,.true.)
   if(a%status/=0)error stop 1
   if(any(shape(a%draws)/=[20,k*d+gdim*k+n*d]))error stop 2
   if(abs(a%draws(1,2)-1.0_dp)>1.0e-12_dp)error stop 3
   if(a%threshold_accept_rate(2)<0.0_dp.or.a%threshold_accept_rate(2)>1.0_dp)error stop 4
   print '(a)','test_ordfactor: PASS'
end program test_ordfactor
