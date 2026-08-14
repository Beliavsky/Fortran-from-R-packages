program test_semiparametric
   use nspmix
   implicit none
   type(nsp_data)::x
   type(nspmix_result)::r
   type(disc_dist)::m
   real(dp)::ni(4),mi(4),ri(4),b0(1)
   integer::grp(6)
   real(dp)::y(6),tr(6),xm(6,1)
   ni=[5.0_dp,6.0_dp,5.0_dp,7.0_dp]; mi=[-1.1_dp,-0.8_dp,2.0_dp,2.2_dp]; ri=[3.0_dp,4.0_dp,3.5_dp,5.0_dp]
   call make_cvps_data(ni,mi,ri,x); call make_disc([-1.0_dp,2.0_dp],d=m); b0=[1.2_dp]
   call cnmms(x,r,b0,m,maxit=2,tol=1.0e-7_dp,ngrid=40,kmax=8)
   if(r%beta(1)<=0.0_dp .or. abs(sum(r%mix%pr)-1.0_dp)>1.0e-10_dp) error stop "cvps"
   grp=[1,1,2,2,3,3]; y=[1.0_dp,2.0_dp,0.0_dp,1.0_dp,3.0_dp,2.0_dp]; tr=4.0_dp
   xm(:,1)=[-1.0_dp,0.5_dp,-0.2_dp,1.0_dp,0.3_dp,1.2_dp]
   call make_mlogit_data(grp,y,tr,xm,x); call make_disc([-1.0_dp,1.0_dp],d=m)
   call cnm_proportions(x,m%pt,r,beta=[0.3_dp],p0=m%pr,maxit=100,tol=1.0e-10_dp)
   if(abs(sum(r%mix%pr)-1.0_dp)>1.0e-10_dp) error stop "mlogit"
   print *, "test_semiparametric: PASS"
end program
