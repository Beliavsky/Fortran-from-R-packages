program test_large_operator
  use gslnls, only : dp, nls_control, nls_result, fit_nls_large_operator, NLS_SUCCESS
  implicit none
  real(dp),parameter::a(4,2)=reshape([1.0_dp,0.0_dp,1.0_dp,2.0_dp, 0.0_dp,1.0_dp,1.0_dp,-1.0_dp],[4,2])
  real(dp),parameter::truth(2)=[2.0_dp,-1.0_dp]
  real(dp)::y(4),start(2)
  type(nls_control)::ctl
  type(nls_result)::fit
  y=matmul(a,truth); start=0.0_dp; ctl=nls_control(); ctl%maxiter=100
  call fit_nls_large_operator(model,jop,y,start,fit,ctl)
  if(fit%status/=NLS_SUCCESS) error stop 'large operator fit failed'
  if(maxval(abs(fit%par-truth))>1.0e-9_dp .or. fit%ssr>1.0e-16_dp) error stop 'large operator mismatch'
  print '(a)', 'PASS test_large_operator'
contains
  subroutine model(par,yhat,ierr)
    real(dp),intent(in)::par(:); real(dp),intent(out)::yhat(:); integer,intent(out)::ierr
    yhat=matmul(a,par); ierr=0
  end subroutine model
  subroutine jop(par,transpose_j,u,v,ierr)
    real(dp),intent(in)::par(:),u(:); logical,intent(in)::transpose_j
    real(dp),intent(out)::v(:); integer,intent(out)::ierr
    if(transpose_j) then; v=matmul(transpose(a),u); else; v=matmul(a,u); end if
    ierr=0
    if(size(par)/=2) ierr=1
  end subroutine jop
end program test_large_operator
