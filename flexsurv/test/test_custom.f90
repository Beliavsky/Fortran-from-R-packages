program test_custom
  use flexsurv
  implicit none
  type(custom_distribution)::d
  type(custom_fit_result)::fit
  type(flexsurv_data)::dat
  real(dp)::t(6),rate
  integer::st(6),fails
  fails=0;t=[0.2_dp,0.5_dp,1.0_dp,1.5_dp,2.0_dp,3.0_dp];st=1
  d%logpdf=>logpdf_exp;d%cdf=>cdf_exp;d%hazard=>haz_exp
  call prepare_survival_data(dat,t,st)
  fit=fit_custom_survival(dat,d,[0.8_dp],positive=[.true.])
  rate=real(size(t),dp)/sum(t)
  if(.not.fit%converged.or.abs(fit%parameters(1)-rate)>3e-5_dp)then
    print *,'custom ',fit%parameters,rate;fails=fails+1
  end if
  if(fails>0)error stop 1
  print *,'test_custom: PASS'
contains
  real(dp) function logpdf_exp(x,p) result(v)
    real(dp),intent(in)::x,p(:);v=log(p(1))-p(1)*x
  end function logpdf_exp
  real(dp) function cdf_exp(x,p) result(v)
    real(dp),intent(in)::x,p(:)
    if(x>1e100_dp)then;v=1.0_dp;else;v=1.0_dp-exp(-p(1)*max(x,0.0_dp));end if
  end function cdf_exp
  real(dp) function haz_exp(x,p) result(v)
    real(dp),intent(in)::x,p(:);v=p(1)+0.0_dp*x
  end function haz_exp
end program test_custom
