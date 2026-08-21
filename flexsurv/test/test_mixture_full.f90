program test_mixture_full
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : flexsurv_data, flexsurv_spec, prepare_survival_data, initialize_spec
  use flexsurv_distributions, only : dist_exponential
  use flexsurv_mixture_full, only : flexsurvmix_full_result, fit_flexsurvmix_full, &
    mix_probability_at, mix_var_louis
  implicit none
  integer, parameter :: n=200
  real(dp) :: time(n),x(n,1),p(2)
  integer :: status(n),ev(n),i
  type(flexsurv_data) :: dat
  type(flexsurv_spec) :: sp(2)
  type(flexsurvmix_full_result) :: fit

  status=1
  do i=1,100
    x(i,1)=-1.0_dp
    if(i<=80)then;ev(i)=1;time(i)=1.0_dp;else;ev(i)=2;time(i)=3.0_dp;end if
  end do
  do i=101,200
    x(i,1)=1.0_dp
    if(i<=120)then;ev(i)=1;time(i)=1.0_dp;else;ev(i)=2;time(i)=3.0_dp;end if
  end do
  call prepare_survival_data(dat,time,status)
  call initialize_spec(sp(1),dist_exponential,n,[1.0_dp])
  call initialize_spec(sp(2),dist_exponential,n,[0.3333333333333333_dp])
  fit=fit_flexsurvmix_full(dat,sp,event=ev,prob_x=x,var_method=mix_var_louis,maxit=100)
  if(.not.fit%converged) error stop 'mixture full did not converge'
  call mix_probability_at(fit,[-1.0_dp],p)
  if(abs(p(2)-0.2_dp)>0.03_dp) error stop 'probability x=-1'
  call mix_probability_at(fit,[1.0_dp],p)
  if(abs(p(2)-0.8_dp)>0.03_dp) error stop 'probability x=1'
  if(.not.allocated(fit%covariance_louis)) error stop 'Louis covariance absent'
  if(size(fit%covariance_louis,1)/=fit%npar) error stop 'Louis covariance shape'
  if(abs(fit%posterior(1,1)-1.0_dp)>1.0e-12_dp) error stop 'known event posterior'
  if(abs(fit%posterior(n,2)-1.0_dp)>1.0e-12_dp) error stop 'known event posterior 2'
  print *, 'test_mixture_full: PASS'
end program test_mixture_full
