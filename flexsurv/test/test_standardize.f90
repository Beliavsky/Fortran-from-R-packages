program test_standardize
  use flexsurv
  implicit none
  type(flexsurv_spec)::sp
  real(dp)::theta(2),s,h,r,exact
  integer::fails
  fails=0
  call initialize_spec(sp,dist_exponential,2,[0.5_dp])
  deallocate(sp%reg(1)%x);allocate(sp%reg(1)%x(2,1));sp%reg(1)%x(:,1)=[0.0_dp,1.0_dp]
  theta=[log(0.5_dp),log(2.0_dp)]
  s=standsurv_survival(sp,theta,1.0_dp);exact=0.5_dp*(exp(-0.5_dp)+exp(-1.0_dp))
  if(abs(s-exact)>1e-12_dp)fails=fails+1
  h=standsurv_hazard(sp,theta,1.0_dp)
  exact=(0.5_dp*exp(-0.5_dp)+1.0_dp*exp(-1.0_dp))/(exp(-0.5_dp)+exp(-1.0_dp))
  if(abs(h-exact)>1e-12_dp)fails=fails+1
  r=standsurv_rmst(sp,theta,2.0_dp)
  exact=0.5_dp*(2.0_dp*(1.0_dp-exp(-1.0_dp))+(1.0_dp-exp(-2.0_dp)))
  if(abs(r-exact)>2e-8_dp)then;print *,'rmst ',r,exact;fails=fails+1;end if
  if(fails>0)error stop 1
  print *,'test_standardize: PASS'
end program test_standardize
