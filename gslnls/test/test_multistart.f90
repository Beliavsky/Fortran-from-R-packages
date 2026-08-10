program test_multistart
  use gslnls, only : dp, nls_control, multistart_result, fit_nls_multistart, NLS_SUCCESS
  implicit none
  real(dp),parameter::xdat(6)=[1.0_dp,2.0_dp,3.0_dp,5.0_dp,7.0_dp,10.0_dp]
  real(dp),parameter::y(6)=[109.0_dp,149.0_dp,149.0_dp,191.0_dp,213.0_dp,224.0_dp]
  real(dp)::ranges(2,2)
  type(nls_control)::ctl
  type(multistart_result)::ans
  ranges(:,1)=[200.0_dp,250.0_dp]; ranges(:,2)=[0.0_dp,1.0_dp]
  ctl=nls_control(); ctl%mstart_n=8; ctl%mstart_q=2; ctl%mstart_s=2
  ctl%mstart_maxstart=30; ctl%mstart_r=1.1_dp; ctl%maxiter=300
  call fit_nls_multistart(model,y,ranges,ans,ctl,jac=jacobian)
  if(ans%fit%status/=NLS_SUCCESS) error stop 'multistart failed'
  if(maxval(abs(ans%fit%par-[213.80940889_dp,0.54723748542_dp]))>5.0e-4_dp) error stop 'multistart mismatch'
  if(abs(ans%fit%ssr-1168.0088766_dp)>1.0e-3_dp) error stop 'multistart rss mismatch'
  if(ans%local_searches<1) error stop 'multistart did no local searches'
  print '(a)', 'PASS test_multistart'
contains
  subroutine model(par,yhat,ierr)
    real(dp),intent(in)::par(:); real(dp),intent(out)::yhat(:); integer,intent(out)::ierr
    yhat=par(1)*(1.0_dp-exp(-par(2)*xdat)); ierr=0
  end subroutine model
  subroutine jacobian(par,j,ierr)
    real(dp),intent(in)::par(:); real(dp),intent(out)::j(:,:); integer,intent(out)::ierr
    real(dp)::e(6); e=exp(-par(2)*xdat); j(:,1)=1.0_dp-e; j(:,2)=par(1)*xdat*e; ierr=0
  end subroutine jacobian
end program test_multistart
