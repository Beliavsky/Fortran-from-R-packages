program test_families
   use nspmix
   implicit none
   type(nsp_data)::x
   type(disc_dist)::m
   real(dp),allocatable::ld(:,:),dt(:,:),db(:,:,:),beta(:)
   real(dp)::ref
   call make_npnorm_data([-1.0_dp,0.0_dp,2.0_dp],x)
   beta=[2.0_dp]; call logd_eval(x,beta,[0.5_dp],ld,dt,db)
   ref=-0.5_dp*(log(2.0_dp*acos(-1.0_dp)*4.0_dp)+((-1.0_dp-0.5_dp)**2)/4.0_dp)
   if(abs(ld(1,1)-ref)>1.0e-13_dp) error stop "normal logd"
   call make_disc([-1.0_dp,2.0_dp],[0.3_dp,0.7_dp],m)
   if(abs(dnpnorm(0.4_dp,m,1.3_dp)-(0.3_dp*exp(-0.5_dp*((0.4_dp+1.0_dp)/1.3_dp)**2) &
      +0.7_dp*exp(-0.5_dp*((0.4_dp-2.0_dp)/1.3_dp)**2))/(1.3_dp*sqrt(2.0_dp*acos(-1.0_dp))))>1.0e-14_dp) error stop "dnpnorm"
   call make_npgeom_data([0.0_dp,1.0_dp,2.0_dp],x)
   call logd_eval(x,[real(dp)::],[0.25_dp],ld,dt,db)
   if(abs(exp(ld(3,1))-0.25_dp*0.75_dp**2)>1.0e-14_dp) error stop "geom"
   call make_npnbinom_data([0.0_dp,3.0_dp],5.0_dp,x)
   call logd_eval(x,[real(dp)::],[0.4_dp],ld,dt,db)
   ref=exp(log_gamma(8.0_dp)-log_gamma(5.0_dp)-log_gamma(4.0_dp)+5.0_dp*log(0.4_dp)+3.0_dp*log(0.6_dp))
   if(abs(exp(ld(2,1))-ref)>1.0e-13_dp) error stop "nbinom"
   print *, "test_families: PASS"
end program
