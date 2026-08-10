program test_algorithms
  use gslnls, only : dp, nls_control, nls_result, fit_nls, NLS_SUCCESS
  use gslnls, only : NLS_LM, NLS_LMACCEL, NLS_DOGLEG, NLS_DDOGLEG, NLS_SUBSPACE2D, NLS_CGST
  implicit none
  real(dp), parameter :: xdat(14) = [77.6_dp,114.9_dp,141.1_dp,190.8_dp,239.9_dp,289.0_dp,332.8_dp, &
       378.4_dp,434.8_dp,477.3_dp,536.8_dp,593.1_dp,689.1_dp,760.0_dp]
  real(dp), parameter :: ydat(14) = [10.07_dp,14.73_dp,17.94_dp,23.93_dp,29.61_dp,35.18_dp,40.02_dp, &
       44.82_dp,50.76_dp,55.05_dp,61.01_dp,66.40_dp,75.47_dp,81.78_dp]
  integer, parameter :: methods(6) = [NLS_LM,NLS_LMACCEL,NLS_DOGLEG,NLS_DDOGLEG,NLS_SUBSPACE2D,NLS_CGST]
  real(dp), parameter :: target(2) = [238.94212918_dp, 5.5015643181e-4_dp]
  real(dp), parameter :: start(2) = [500.0_dp,1.0e-4_dp]
  type(nls_control) :: ctl
  type(nls_result) :: fit
  integer :: k

  do k=1,size(methods)
    ctl=nls_control(); ctl%algorithm=methods(k); ctl%maxiter=500
    call fit_nls(model,ydat,start,fit,ctl,jac=jacobian)
    if(fit%status/=NLS_SUCCESS) error stop 'algorithm did not converge'
    if(maxval(abs(fit%par-target))>5.0e-4_dp) error stop 'algorithm parameter mismatch'
    if(abs(fit%ssr-0.1245514_dp)>2.0e-5_dp) error stop 'algorithm rss mismatch'
  end do
  print '(a)', 'PASS test_algorithms'
contains
  subroutine model(par,yhat,ierr)
    real(dp),intent(in)::par(:)
    real(dp),intent(out)::yhat(:)
    integer,intent(out)::ierr
    yhat=par(1)*(1.0_dp-exp(-par(2)*xdat)); ierr=0
  end subroutine model
  subroutine jacobian(par,j,ierr)
    real(dp),intent(in)::par(:)
    real(dp),intent(out)::j(:,:)
    integer,intent(out)::ierr
    real(dp)::e(size(xdat))
    e=exp(-par(2)*xdat)
    j(:,1)=1.0_dp-e
    j(:,2)=par(1)*xdat*e
    ierr=0
  end subroutine jacobian
end program test_algorithms
