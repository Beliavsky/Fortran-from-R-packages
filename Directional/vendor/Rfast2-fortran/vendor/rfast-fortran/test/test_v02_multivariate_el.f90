program test_v02_multivariate_el
   use rfast
   implicit none
   real(dp)::x(6,2),xl(6,2),dir(6,3),invx(6,2),v(7),x1(5),x2(5),mu2(2)
   type(multivariate_mle_result)::m
   type(el_result)::e
   type(el_two_sample_result)::e2
   type(mv_el_result)::me
   integer::fail
   fail=0
   x=reshape([ -2.0_dp,-1.0_dp, 0.0_dp,1.0_dp,2.0_dp,0.0_dp, &
               1.0_dp, 0.0_dp,-1.0_dp,0.0_dp,1.0_dp,-1.0_dp ],shape(x))
   m=mvnorm_mle(x)
   if(maxval(abs(m%location-colmeans(x)))>1e-12_dp)fail=fail+1
   xl=exp(0.2_dp*x+1.0_dp);m=mvlognormal_mle(xl);if(m%status/=0.or.any(m%mean_original<=0.0_dp))fail=fail+1
   m=mvt_mle(x,5.0_dp);if(m%status/=0.or..not.allocated(m%scatter))fail=fail+1
   dir=reshape([0.2_dp,0.3_dp,0.5_dp, 0.25_dp,0.25_dp,0.5_dp, 0.3_dp,0.2_dp,0.5_dp, &
                0.15_dp,0.35_dp,0.5_dp, 0.22_dp,0.28_dp,0.5_dp, 0.28_dp,0.32_dp,0.4_dp],shape(dir),order=[2,1])
   m=dirichlet_mle(dir);if(m%status/=0.or.any(m%param<=0.0_dp))fail=fail+1
   invx=reshape([0.5_dp,1.0_dp,1.5_dp,0.8_dp,1.2_dp,2.0_dp, 1.0_dp,0.7_dp,1.2_dp,1.8_dp,0.9_dp,1.4_dp],shape(invx))
   m=inverse_dirichlet_mle(invx);if(m%status/=0.or.any(m%param<=0.0_dp))fail=fail+1
   v=[-3.0_dp,-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp,3.0_dp];e=el_test1(v,0.0_dp)
   if(abs(e%lambda)>1e-9_dp.or.abs(e%statistic)>1e-9_dp.or.maxval(abs(e%weights-1.0_dp/7.0_dp))>1e-9_dp)fail=fail+1
   e=eel_test1(v,0.0_dp);if(abs(e%lambda)>1e-9_dp.or.abs(e%statistic)>1e-9_dp)fail=fail+1
   x1=[-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp];x2=x1;e2=el_test2(x1,x2);if(e2%statistic>1e-8_dp)fail=fail+1
   mu2=0.0_dp;me=mv_eeltest1(x,mu2);if(me%status/=0.or.abs(me%statistic)>1e-7_dp)fail=fail+1
   me=mv_eeltest2(x,x);if(me%status/=0.or.abs(me%statistic)>1e-7_dp)fail=fail+1
   if(fail==0)then
      print '(a)','test_v02_multivariate_el: PASS'
   else
      print '(a,i0)','test_v02_multivariate_el: FAIL ',fail
      error stop 1
   end if
end program test_v02_multivariate_el
