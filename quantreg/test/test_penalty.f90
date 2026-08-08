program test_penalty
  use quantreg, only : dp, rq_result, rq_fit_lasso, rq_fit_scad
  implicit none
  real(dp) :: x(9,2), y(9), lam(2)
  type(rq_result) :: lfit, sfit
  integer :: i
  do i=1,9
    x(i,1)=1.0_dp
    x(i,2)=real(i-5,dp)
  end do
  y = 2.0_dp + 3.0_dp*x(:,2)
  y(1)=y(1)-4.0_dp
  y(9)=y(9)+4.0_dp
  lam=[0.0_dp,0.5_dp]
  call rq_fit_lasso(x,y,0.5_dp,lam,lfit)
  if (lfit%info /= 0) error stop 'lasso info'
  if (abs(lfit%coefficients(1)-2.0_dp) > 1.0e-4_dp) error stop 'lasso intercept'
  call rq_fit_scad(x,y,0.5_dp,lam,sfit)
  if (sfit%info /= 0) error stop 'scad info'
  if (abs(sfit%coefficients(1)-2.0_dp) > 1.0e-3_dp) error stop 'scad intercept'
  print *, 'test_penalty: PASS'
end program
