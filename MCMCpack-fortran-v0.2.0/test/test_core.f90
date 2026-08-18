program test_core
   use mcmcpack
   implicit none
   real(dp) :: x3(3),a3(3),d,mat(3,3),back(3,3),w1(1,1),s1(1,1),p0
   real(dp) :: draws(6000), md(3000,3), meanv
   integer :: status,vals_x
   call set_seed(12345)
   if(abs(dinvgamma(1.0_dp,2.0_dp,3.0_dp)-9.0_dp*exp(-3.0_dp))>1.0e-12_dp) error stop 1
   x3=[0.2_dp,0.3_dp,0.5_dp];a3=1.0_dp
   if(abs(ddirichlet(x3,a3)-2.0_dp)>1.0e-12_dp) error stop 2
   p0=dnoncenhypergeom(2,5,7,4,1.0_dp)
   if(abs(p0-(10.0_dp*21.0_dp/495.0_dp))>1.0e-12_dp) error stop 3
   mat=reshape([1.0_dp,2.0_dp,3.0_dp,2.0_dp,4.0_dp,5.0_dp,3.0_dp,5.0_dp,7.0_dp],[3,3])
   back=xpnd(vech(mat),3)
   if(maxval(abs(back-mat))>1.0e-12_dp) error stop 4
   w1(1,1)=2.0_dp;s1(1,1)=1.0_dp;d=dwish(w1,3.0_dp,s1)
   if(.not.(d>0.0_dp)) error stop 5
   call mc_binomial_beta(3,12,1.0_dp,1.0_dp,draws,status)
   if(status/=0) error stop 6
   meanv=sum(draws)/real(size(draws),dp)
   if(abs(meanv-4.0_dp/14.0_dp)>0.02_dp) error stop 7
   call mc_multinom_dirichlet([10,20,30],[1.0_dp,1.0_dp,1.0_dp],md,status)
   if(status/=0.or.maxval(abs(sum(md,dim=2)-1.0_dp))>1.0e-12_dp) error stop 8
   vals_x=rnoncenhypergeom(5,7,4,2.0_dp)
   if(vals_x<0.or.vals_x>4) error stop 9
   print *, 'test_core: PASS'
end program test_core
