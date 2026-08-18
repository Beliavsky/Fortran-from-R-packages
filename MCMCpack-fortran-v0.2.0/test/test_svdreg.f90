program test_svdreg
   use mcmcpack
   implicit none
   real(dp)::y(3),a(2,3),d(3,3),f(3,3),tau(3),g0(3),c0(3),d0(3),w0(3)
   logical::miss(3)
   type(mcmc_result)::r
   call set_seed(86420)
   y=[1.0_dp,0.0_dp,-0.5_dp];miss=[.false.,.true.,.false.]
   a=reshape([1.0_dp,0.0_dp,0.5_dp, 0.0_dp,1.0_dp,0.5_dp],[2,3],order=[2,1])
   a(1,:)=[1.0_dp,0.0_dp,0.5_dp];a(2,:)=[0.0_dp,1.0_dp,0.5_dp]
   d=0.0_dp;d(1,1)=1.0_dp;d(2,2)=1.5_dp;d(3,3)=2.0_dp
   f=0.0_dp;f(1,1)=1.0_dp;f(2,2)=1.0_dp;f(3,3)=1.0_dp
   tau=1.0_dp;g0=0.0_dp;c0=2.0_dp;d0=2.0_dp;w0=0.3_dp
   r=mcmc_svdreg(y,miss,a,d,f,tau,g0,2.0_dp,2.0_dp,c0,d0,w0,10,20,2,.true.)
   if(r%status/=0) error stop 'svdreg status'
   if(size(r%draws,1)/=10.or.size(r%draws,2)/=1+3+3+1+2) error stop 'svdreg dimensions'
   if(any(r%draws(:,1+3+1:1+3+3)<=0.0_dp)) error stop 'svdreg tau'
   print '(a)','test_svdreg: PASS'
end program
