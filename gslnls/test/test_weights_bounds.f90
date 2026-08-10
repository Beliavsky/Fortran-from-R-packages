program test_weights_bounds
  use gslnls, only : dp, nls_control, nls_result, fit_nls, NLS_SUCCESS
  implicit none
  real(dp), parameter :: xdat(14) = [77.6_dp,114.9_dp,141.1_dp,190.8_dp,239.9_dp,289.0_dp,332.8_dp, &
       378.4_dp,434.8_dp,477.3_dp,536.8_dp,593.1_dp,689.1_dp,760.0_dp]
  real(dp), parameter :: ydat(14) = [10.07_dp,14.73_dp,17.94_dp,23.93_dp,29.61_dp,35.18_dp,40.02_dp, &
       44.82_dp,50.76_dp,55.05_dp,61.01_dp,66.40_dp,75.47_dp,81.78_dp]
  real(dp), parameter :: start(2)=[500.0_dp,1.0e-4_dp]
  real(dp) :: w(14), wm(14,14), lo(2), hi(2)
  type(nls_result)::a,b,c
  type(nls_control)::ctl
  integer::i
  w=100.0_dp; wm=0.0_dp
  do i=1,14; wm(i,i)=100.0_dp; end do
  ctl=nls_control(); ctl%maxiter=500
  call fit_nls(model,ydat,start,a,ctl,jac=jacobian,weights=w)
  call fit_nls(model,ydat,start,b,ctl,jac=jacobian,weight_matrix=wm)
  if(a%status/=NLS_SUCCESS .or. b%status/=NLS_SUCCESS) error stop 'weighted fit failed'
  if(maxval(abs(a%par-b%par))>1.0e-6_dp) error stop 'vector/matrix weights disagree'
  if(abs(a%ssr-12.4551389_dp)>2.0e-4_dp) error stop 'weighted rss mismatch'
  lo=[250.0_dp,0.0_dp]; hi=[huge(1.0_dp),1.0_dp]
  call fit_nls(model,ydat,start,c,ctl,jac=jacobian,lower=lo,upper=hi)
  if(c%status/=NLS_SUCCESS) error stop 'bounded fit failed'
  if(abs(c%par(1)-250.0_dp)>1.0e-8_dp .or. abs(c%par(2)-5.22025681e-4_dp)>5.0e-7_dp) &
       error stop 'bounded result mismatch'
  print '(a)', 'PASS test_weights_bounds'
contains
  subroutine model(par,yhat,ierr)
    real(dp),intent(in)::par(:); real(dp),intent(out)::yhat(:); integer,intent(out)::ierr
    yhat=par(1)*(1.0_dp-exp(-par(2)*xdat)); ierr=0
  end subroutine model
  subroutine jacobian(par,j,ierr)
    real(dp),intent(in)::par(:); real(dp),intent(out)::j(:,:); integer,intent(out)::ierr
    real(dp)::e(14); e=exp(-par(2)*xdat); j(:,1)=1.0_dp-e; j(:,2)=par(1)*xdat*e; ierr=0
  end subroutine jacobian
end program test_weights_bounds
