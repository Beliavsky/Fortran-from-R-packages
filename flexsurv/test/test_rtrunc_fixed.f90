program test_rtrunc_fixed
  use flexsurv_kinds, only : dp
  use flexsurv_distributions, only : dist_gamma
  use flexsurv_rtrunc, only : flexsurvrtrunc_result, fit_flexsurvrtrunc, rtrunc_final, rtrunc_parameter_draws
  implicit none
  type(flexsurvrtrunc_result)::fit
  real(dp)::t(8),ti(8),rt(8)
  real(dp),allocatable::draws(:,:)
  logical::fixp(2)
  integer::fails
  fails=0;t=[0.5_dp,0.8_dp,1.1_dp,1.4_dp,1.8_dp,2.2_dp,2.6_dp,3.0_dp];ti=0.0_dp;rt=10.0_dp
  fixp=[.true.,.false.]
  fit=fit_flexsurvrtrunc(t,ti,rt,10.0_dp,dist_gamma,[2.0_dp,1.0_dp],0.2_dp, &
    method=rtrunc_final,fixed_theta=.true.,fixed_params=fixp,maxit=300)
  if(.not.fit%converged)then;print *,'FAIL convergence';fails=fails+1;end if
  if(abs(fit%parameters(1)-2.0_dp)>1.0e-12_dp)then;print *,'FAIL fixed shape';fails=fails+1;end if
  if(abs(fit%covariance(1,1))>1.0e-15_dp.or.abs(fit%covariance(3,3))>1.0e-15_dp)then
    print *,'FAIL fixed covariance';fails=fails+1
  end if
  call rtrunc_parameter_draws(fit,20,draws,seed=44)
  if(maxval(abs(draws(1,:)-2.0_dp))>1.0e-12_dp)then;print *,'FAIL fixed draws';fails=fails+1;end if
  if(fails==0)then;print *,'test_rtrunc_fixed: PASS';else;print *,'test_rtrunc_fixed: FAIL',fails;error stop 1;end if
end program test_rtrunc_fixed
