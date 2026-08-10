program robust_line
  use gslnls, only : dp, nls_control, nls_result, fit_nls, default_loss, LOSS_HUBER
  implicit none
  real(dp)::x(21),y(21)
  type(nls_control)::ctl
  type(nls_result)::fit
  integer::i
  do i=1,21; x(i)=0.1_dp*real(i-1,dp); end do
  y=1.0_dp+2.0_dp*x; y(11)=20.0_dp
  ctl=nls_control()
  call fit_nls(model,y,[0.0_dp,0.0_dp],fit,ctl,jac=jacobian,loss=default_loss(LOSS_HUBER))
  print '(a,2f14.8)', 'Huber line parameters:',fit%par
  print '(a,i0)', 'IRLS iterations: ',fit%irls_iterations
contains
  subroutine model(par,yhat,ierr)
    real(dp),intent(in)::par(:); real(dp),intent(out)::yhat(:); integer,intent(out)::ierr
    yhat=par(1)+par(2)*x; ierr=0
  end subroutine model
  subroutine jacobian(par,j,ierr)
    real(dp),intent(in)::par(:); real(dp),intent(out)::j(:,:); integer,intent(out)::ierr
    j(:,1)=1.0_dp; j(:,2)=x; ierr=0
    if(size(par)/=2) ierr=1
  end subroutine jacobian
end program robust_line
