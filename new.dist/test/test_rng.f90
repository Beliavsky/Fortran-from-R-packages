program test_rng
   use new_dist
   implicit none
   integer,parameter::n=1000
   real(dp)::xr(n)
   integer::xi(n),fails
   fails=0
   call set_new_dist_seed(12345)
   call rEPd(xr,2.0_dp,3.0_dp); if(any(xr<=0.0_dp)) fails=fails+1
   call rLd(xr,2.0_dp); if(any(xr<0.0_dp)) fails=fails+1
   call rRA(xr,1.0_dp); if(any(xr<0.0_dp)) fails=fails+1
   call rbwd(xr,2.0_dp,1.0_dp,2.0_dp); if(any(xr<0.0_dp)) fails=fails+1
   call rdLd1(xi,1.0_dp); if(any(xi<0)) fails=fails+1
   call rdLd2(xi,1.0_dp); if(any(xi<0)) fails=fails+1
   call rgld(xr,2.0_dp,3.0_dp,4.0_dp); if(any(xr<0.0_dp)) fails=fails+1
   call rkd(xr,2.0_dp,3.0_dp); if(any(xr<0.0_dp).or.any(xr>1.0_dp)) fails=fails+1
   call rmd(xr,2.0_dp); if(any(xr<0.0_dp)) fails=fails+1
   call romd(xr,0.2_dp); if(any(xr<0.0_dp)) fails=fails+1
   call rpldd(xr,2.0_dp,3.0_dp,4.0_dp)
   call rsgrd(xr,2.0_dp,1.0_dp,3.0_dp); if(any(xr<0.0_dp)) fails=fails+1
   call rsod(xr,1.0_dp,2.0_dp); if(any(xr<0.0_dp).or.any(xr>1.0_dp)) fails=fails+1
   call rtpmd(xr,1.0_dp,2.0_dp); if(any(xr<0.0_dp)) fails=fails+1
   call rtprd(xr,1.0_dp,1.0_dp); if(any(xr<1.0_dp)) fails=fails+1
   call rugd(xi,0.5_dp); if(any(xi<1)) fails=fails+1
   call ruigd(xr,1.0_dp,1.0_dp); if(any(xr<=0.0_dp)) fails=fails+1
   call rwgd(xi,0.2_dp,3.0_dp); if(any(xi<1)) fails=fails+1
   if(fails/=0) then; print *,'test_rng: FAIL',fails; error stop 1; end if
   print *,'test_rng: PASS'
end program test_rng
