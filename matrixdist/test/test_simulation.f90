program test_simulation
   use r_compat, only: dp
   use matrixdist_simulation
   use matrixdist_transformations
   implicit none
   integer,parameter::n=30000
   real(dp)::a(1),s(1,1),m
   real(dp),allocatable::x(:),r(:,:),bm(:,:)
   integer,allocatable::k(:)
   a=1.0_dp
   s=-2.0_dp
   x=rphasetype(n,a,s)
   m=sum(x)/n
   if(abs(m-0.5_dp)>0.02_dp) then
   print *,'ph sim',m
   error stop 1
   end if
   s=0.6_dp
   k=rdphasetype(n,a,s)
   m=real(sum(k),dp)/n
   if(abs(m-2.5_dp)>0.08_dp) then
   print *,'dph sim',m
   error stop 1
   end if
   s=-2.0_dp
   biv_sim: block
      real(dp)::s11(1,1),s12(1,1),s22(1,1)
      s11=-2.0_dp
      s12=2.0_dp
      s22=-3.0_dp
      bm=rbivph(n,a,s11,s12,s22)
      if(abs(sum(bm(:,1))/n-0.5_dp)>0.02_dp) error stop 'bivph sim mean1'
      if(abs(sum(bm(:,2))/n-1.0_dp/3.0_dp)>0.02_dp) error stop 'bivph sim mean2'
   end block biv_sim
   r=random_reward(5,3)
   if(maxval(abs(sum(r,dim=2)-1.0_dp))>1e-12_dp) error stop 'reward rows'
   print *, 'test_simulation: PASS'
end program
