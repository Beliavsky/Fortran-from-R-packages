program test_spline_fit
  use flexsurv
  implicit none
  type(flexsurv_data)::dat
  type(flexsurvspline_result)::fit
  real(dp)::t(12),kn(2)
  integer::st(12),i,fails
  fails=0
  do i=1,12;t(i)=-log(1.0_dp-(real(i,dp)-0.5_dp)/12.0_dp);end do
  st=1;call prepare_survival_data(dat,t,st)
  kn=[log(0.05_dp),log(20.0_dp)]
  fit=fit_flexsurvspline(dat,kn,spline_scale_hazard,spline_time_log, &
    gamma_init=[0.0_dp,1.0_dp],maxit=300,tol=1e-7_dp)
  if(.not.fit%converged)then;print *,'spline fit status ',fit%status;fails=fails+1;end if
  if(abs(fit%model%gamma(1))>0.15_dp.or.abs(fit%model%gamma(2)-1.0_dp)>0.15_dp)then
    print *,'gamma ',fit%model%gamma;fails=fails+1
  end if
  if(fails>0)error stop 1
  print *,'test_spline_fit: PASS'
end program test_spline_fit
