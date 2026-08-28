program test_spline_interactions
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : flexsurv_data, parameter_regression, prepare_survival_data
  use flexsurv_spline, only : survspline_model, survspline_quantile, spline_scale_hazard, spline_time_log
  use flexsurv_spline_interactions, only : flexsurvspline_interaction_result, &
    fit_flexsurvspline_interactions, spline_interaction_gamma
  implicit none
  integer,parameter::n=160
  type(flexsurv_data)::data
  type(parameter_regression)::reg(2)
  type(flexsurvspline_interaction_result)::fit
  type(survspline_model)::m
  real(dp)::time(n),x(n),knots(2),gtrue(2),p,beta_true
  real(dp),allocatable::ghat(:)
  integer::i,fails,status(n)

  fails=0
  knots=[-3.0_dp,3.0_dp]
  gtrue=[-1.0_dp,1.25_dp]
  beta_true=0.45_dp
  allocate(m%knots(2),m%gamma(2))
  m%knots=knots
  m%scale=spline_scale_hazard
  m%timescale=spline_time_log
  do i=1,n
    if(i<=n/2)then
    x(i)=-0.5_dp
    else
    x(i)=0.5_dp
    end if
    p=(real(mod(i-1,n/2),dp)+0.5_dp)/real(n/2,dp)
    m%gamma=gtrue
    m%gamma(2)=m%gamma(2)+beta_true*x(i)
    time(i)=survspline_quantile(m,p)
  end do
  status=1
  call prepare_survival_data(data,time,status=status)
  allocate(reg(1)%x(n,0),reg(2)%x(n,1))
  reg(2)%x(:,1)=x
  fit=fit_flexsurvspline_interactions(data,knots,spline_scale_hazard,spline_time_log,reg, &
    gamma_init=[-0.9_dp,1.2_dp],beta_init=[0.3_dp],maxit=500)
  if(.not.fit%converged)then
  print *,'FAIL convergence',fit%status,fit%iterations,maxval(abs(fit%gradient))
  fails=fails+1
  end if
  if(abs(fit%model%gamma(1)-gtrue(1))>0.18_dp)then
  print *,'FAIL gamma0',fit%model%gamma(1)
  fails=fails+1
  end if
  if(abs(fit%model%gamma(2)-gtrue(2))>0.18_dp)then
  print *,'FAIL gamma1',fit%model%gamma(2)
  fails=fails+1
  end if
  if(abs(fit%effect(2)%beta(1)-beta_true)>0.20_dp)then
    print *,'FAIL gamma1 covariate',fit%effect(2)%beta(1)
    fails=fails+1
  end if
  ghat=spline_interaction_gamma(reg,fit%theta,n,2)
  if(abs(ghat(2)-(fit%model%gamma(2)+0.5_dp*fit%effect(2)%beta(1)))>1.0e-12_dp)then
    print *,'FAIL row gamma'
    fails=fails+1
  end if
  if(.not.allocated(fit%covariance))then
  print *,'FAIL covariance'
  fails=fails+1
  end if
  if(fails==0)then
    print *,'test_spline_interactions: PASS'
  else
    print *,'test_spline_interactions: FAIL',fails
    error stop 1
  end if
end program test_spline_interactions
