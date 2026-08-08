program constrained_regression
  use quantreg, only : dp, rq_result, rq_fit_fnc
  implicit none
  real(dp) :: x(5,2), y(5), rmat(1,2), rhs(1)
  type(rq_result) :: fit

  x(:,1)=1.0_dp
  x(:,2)=[0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp]
  y=[1.0_dp,2.0_dp,4.0_dp,5.0_dp,7.0_dp]
  rmat=0.0_dp
  rmat(1,2)=1.0_dp
  rhs=2.0_dp
  call rq_fit_fnc(x,y,rmat,rhs,0.5_dp,fit)
  print '(a,2f14.7)', 'constrained coefficients: ', fit%coefficients
end program constrained_regression
