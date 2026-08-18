program test_hglm
   use mcmcpack
   implicit none
   integer,parameter::n=12
   real(dp)::x(n,2),w(n,1),b0(2),bs(2),bst(3,1),vb(1,1),vbeta(2,2),rmat(1,1),z0(n)
   integer::grp(n),yb(n),yp(n),i
   type(hglm_result)::a,b
   x(:,1)=1.0_dp;do i=1,n;x(i,2)=real(i-6,dp)/5.0_dp;grp(i)=1+mod(i-1,3);end do
   w=1.0_dp;yb=merge(1,0,x(:,2)>0.0_dp);yp=[0,1,0,1,2,1,2,3,2,3,4,3]
   b0=0.0_dp;bs=0.0_dp;bst=0.0_dp;vb=1.0_dp;vbeta=0.0_dp;vbeta(1,1)=10.0_dp;vbeta(2,2)=10.0_dp;rmat=1.0_dp
   z0=merge(0.7_dp,-0.7_dp,yb==1)
   call set_seed(811)
   a=mcmc_hlogit(yb,x,w,grp,bs,bst,vb,1.0_dp,z0,b0,vbeta,3.0_dp,rmat,2.0_dp,2.0_dp,20,40,2)
   if(a%status/=0.or.any(shape(a%draws)/=[20,8]))error stop 1
   if(any(a%prediction<=0.0_dp).or.any(a%prediction>=1.0_dp))error stop 2
   z0=log(real(yp,dp)+0.5_dp)
   call set_seed(812)
   b=mcmc_hpoisson(yp,x,w,grp,bs,bst,vb,1.0_dp,z0,b0,vbeta,3.0_dp,rmat,2.0_dp,2.0_dp,20,40,2)
   if(b%status/=0.or.any(shape(b%draws)/=[20,8]))error stop 3
   if(any(b%prediction<=0.0_dp))error stop 4
   print '(a)','test_hglm: PASS'
end program test_hglm
