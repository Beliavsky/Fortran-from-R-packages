program demo_flexsurv
  use flexsurv
  implicit none
  type(flexsurv_data)::dat
  type(flexsurv_spec)::sp
  type(flexsurv_result)::fit
  real(dp)::t(10)
  integer::st(10)
  t=[0.2_dp,0.4_dp,0.7_dp,0.9_dp,1.1_dp,1.4_dp,1.8_dp,2.0_dp,2.4_dp,3.0_dp]
  st=1
  call prepare_survival_data(dat,t,st)
  call initialize_spec(sp,dist_weibull,size(t),[1.0_dp,1.0_dp])
  fit=fit_flexsurvreg(dat,sp)
  print '(a,l1)','converged: ',fit%converged
  print '(a,*(f10.5,1x))','weibull shape/scale: ',fit%base
  print '(a,f12.6)','log-likelihood: ',fit%loglik
  print '(a,f10.6)','S(1): ',predict_survival(sp,fit%theta,1,1.0_dp)
end program demo_flexsurv
