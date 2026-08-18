program test_weighted
  use quantreg, only : dp, rq_result, rq_wfit_fnb, check_loss_sum
  implicit none
  real(dp) :: x(7,2), y(7), w(7), r(3)
  type(rq_result) :: fit
  integer :: i

  do i=1,7
    x(i,1)=1.0_dp
    x(i,2)=real(i-4,dp)
  end do
  y=1.25_dp-0.75_dp*x(:,2)
  w=[1.0_dp,2.0_dp,0.5_dp,3.0_dp,1.0_dp,1.5_dp,2.5_dp]
  call rq_wfit_fnb(x,y,w,0.5_dp,fit)
  if (fit%info/=0) error stop 'weighted info'
  if (maxval(abs(fit%coefficients-[1.25_dp,-0.75_dp]))>1.0e-6_dp) error stop 'weighted coef'
  r=[-2.0_dp,0.0_dp,3.0_dp]
  if (abs(check_loss_sum(r,0.25_dp)-2.25_dp)>1.0e-12_dp) error stop 'check loss'
  print *, 'test_weighted: PASS'
end program test_weighted
