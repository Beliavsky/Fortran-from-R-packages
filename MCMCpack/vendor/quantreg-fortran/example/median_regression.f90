program median_regression
  use quantreg, only : dp, rq_result, rq_fit_fnb
  implicit none
  real(dp) :: x(5,2), y(5)
  type(rq_result) :: fit

  x(:,1)=1.0_dp
  x(:,2)=[0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp]
  y=[1.0_dp,2.0_dp,4.0_dp,5.0_dp,7.0_dp]
  call rq_fit_fnb(x,y,0.5_dp,fit)
  print '(a,2f14.7)', 'coefficients: ', fit%coefficients
  print '(a,i0)', 'iterations: ', fit%iterations
end program median_regression
