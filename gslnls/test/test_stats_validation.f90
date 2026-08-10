program test_stats_validation
  use gslnls, only : dp, nls_control, nls_result, fit_nls, nls_hatvalues, nls_cooks_distance, nls_loglik
  use gslnls, only : NLS_SUCCESS, NLS_BAD_INPUT
  implicit none
  real(dp),parameter::xdat(5)=[0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp]
  real(dp),parameter::y(5)=[1.0_dp,3.0_dp,5.0_dp,7.0_dp,9.0_dp]
  type(nls_result)::fit,bad
  type(nls_control)::ctl
  real(dp)::h(5),cook(5),lo(2),hi(2),ll
  integer::rank,ierr
  ctl=nls_control(); ctl%store_trace=.true.
  call fit_nls(model,y,[0.0_dp,0.0_dp],fit,ctl,jac=jacobian)
  if(fit%status/=NLS_SUCCESS) error stop 'linear stats fit failed'
  call nls_hatvalues(fit%jacobian,h,rank)
  if(rank/=2 .or. abs(sum(h)-2.0_dp)>1.0e-10_dp) error stop 'hat values wrong'
  call nls_cooks_distance(fit%residual,fit%jacobian,max(fit%sigma,1.0e-12_dp),cook,ierr)
  if(ierr/=0 .or. any(cook<0.0_dp)) error stop 'Cook distance failed'
  ll=nls_loglik(max(fit%ssr,1.0e-20_dp),5)
  if(ll<=-huge(1.0_dp)/2.0_dp) error stop 'loglik failed'
  if(.not.allocated(fit%par_trace) .or. .not.allocated(fit%ssr_trace)) error stop 'trace missing'
  lo=[2.0_dp,0.0_dp]; hi=[1.0_dp,5.0_dp]
  call fit_nls(model,y,[0.0_dp,0.0_dp],bad,ctl,lower=lo,upper=hi)
  if(bad%status/=NLS_BAD_INPUT) error stop 'invalid bounds not rejected'
  print '(a)', 'PASS test_stats_validation'
contains
  subroutine model(par,yhat,ierr)
    real(dp),intent(in)::par(:); real(dp),intent(out)::yhat(:); integer,intent(out)::ierr
    yhat=par(1)+par(2)*xdat; ierr=0
  end subroutine model
  subroutine jacobian(par,j,ierr)
    real(dp),intent(in)::par(:); real(dp),intent(out)::j(:,:); integer,intent(out)::ierr
    j(:,1)=1.0_dp; j(:,2)=xdat; ierr=0
    if(size(par)/=2) ierr=1
  end subroutine jacobian
end program test_stats_validation
