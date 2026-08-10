program test_robust
  use gslnls, only : dp, nls_control, nls_result, nls_loss, default_loss, fit_nls, LOSS_HUBER, NLS_SUCCESS
  implicit none
  real(dp)::xdat(21),y(21),start(2),eols,erob
  type(nls_control)::ctl
  type(nls_result)::ols,rob
  type(nls_loss)::loss
  integer::i
  do i=1,21; xdat(i)=0.1_dp*real(i-1,dp); end do
  y=1.0_dp+2.0_dp*xdat; y(11)=20.0_dp; start=[0.0_dp,0.0_dp]
  ctl=nls_control(); ctl%maxiter=200
  call fit_nls(model,y,start,ols,ctl,jac=jacobian)
  loss=default_loss(LOSS_HUBER)
  call fit_nls(model,y,start,rob,ctl,jac=jacobian,loss=loss)
  if(ols%status/=NLS_SUCCESS .or. rob%status/=NLS_SUCCESS) error stop 'robust fit failed'
  eols=maxval(abs(ols%par-[1.0_dp,2.0_dp])); erob=maxval(abs(rob%par-[1.0_dp,2.0_dp]))
  if(erob>=1.0e-3_dp .or. erob>=eols/100.0_dp) error stop 'Huber IRLS did not reject outlier'
  if(rob%irls_iterations<1 .or. .not.allocated(rob%irls_weights)) error stop 'IRLS bookkeeping missing'
  print '(a)', 'PASS test_robust'
contains
  subroutine model(par,yhat,ierr)
    real(dp),intent(in)::par(:); real(dp),intent(out)::yhat(:); integer,intent(out)::ierr
    yhat=par(1)+par(2)*xdat; ierr=0
  end subroutine model
  subroutine jacobian(par,j,ierr)
    real(dp),intent(in)::par(:); real(dp),intent(out)::j(:,:); integer,intent(out)::ierr
    j(:,1)=1.0_dp; j(:,2)=xdat; ierr=0
    if(size(par)<1) ierr=1
  end subroutine jacobian
end program test_robust
