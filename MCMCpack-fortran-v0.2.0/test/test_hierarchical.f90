program test_hierarchical
   use mcmcpack
   implicit none
   type(hregress_result) :: r
   integer,parameter :: n=12,p=2,q=1,ng=3
   real(dp)::y(n),x(n,p),w(n,q),bst(p),bb(ng,q),vb(q,q),mub(p),vbeta(p,p),rmat(q,q)
   integer::grp(n),i
   call set_seed(24680)
   do i=1,n
      grp(i)=1+mod(i-1,ng);x(i,1)=1.0_dp;x(i,2)=real(i-6,dp)/5.0_dp;w(i,1)=1.0_dp
      y(i)=1.0_dp+0.7_dp*x(i,2)+0.2_dp*real(grp(i)-2,dp)+0.05_dp*real(mod(i,2),dp)
   end do
   bst=0.0_dp;bb=0.0_dp;vb=1.0_dp;mub=0.0_dp;vbeta=0.0_dp;vbeta(1,1)=10.0_dp;vbeta(2,2)=10.0_dp;rmat=1.0_dp
   r=mcmc_hregress(y,x,w,grp,bst,bb,vb,1.0_dp,mub,vbeta,3.0_dp,rmat,2.0_dp,1.0_dp,20,40,2)
   if(r%status/=0) error stop 'hregress status'
   if(size(r%draws,1)/=20 .or. size(r%draws,2)/=8) error stop 'hregress dimensions'
   if(any(.not.(r%draws==r%draws))) error stop 'hregress NaN'
   if(size(r%y_pred)/=n) error stop 'hregress prediction'
   print '(a)', 'test_hierarchical: PASS'
end program
