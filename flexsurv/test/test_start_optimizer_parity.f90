program test_start_optimizer_parity
  use flexsurv_kinds, only : dp
  use flexsurv_fit
  use flexsurv_distributions, only : dist_exponential
  implicit none
  type(flexsurv_data)::dat
  type(flexsurv_spec)::sp
  type(flexsurv_result)::fb,fn
  real(dp),allocatable::z0(:)
  real(dp)::t(8),truth
  integer::st(8)
  t=[0.2_dp,0.4_dp,0.7_dp,1.0_dp,1.3_dp,1.7_dp,2.1_dp,2.9_dp];st=1
  truth=real(size(t),dp)/sum(t)
  call prepare_survival_data(dat,t,st)
  call initialize_spec(sp,dist_exponential,size(t),[1.0_dp])
  z0=source_initial_theta(dat,sp)
  if(abs(exp(z0(1))-truth)>2.0e-5_dp)error stop 'source exp start'
  fb=fit_flexsurvreg(dat,sp,source_start=.true.,optim_method=fs_optim_bfgs)
  fn=fit_flexsurvreg(dat,sp,source_start=.true.,optim_method=fs_optim_nelder_mead)
  if(.not.fb%converged.or..not.fn%converged)error stop 'optimizer convergence'
  if(abs(fb%base(1)-truth)>2.0e-5_dp)error stop 'bfgs exp mle'
  if(abs(fn%base(1)-truth)>2.0e-4_dp)error stop 'simplex exp mle'
  if(abs(fb%loglik-fn%loglik)>1.0e-6_dp)error stop 'optimizer agreement'
  print *,'test_start_optimizer_parity: PASS'
end program
