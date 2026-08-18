program test_nonlinear
  use quantreg, only : dp, nlrq_result, nlrq_fit
  implicit none
  real(dp) :: y(9), theta0(2)
  type(nlrq_result) :: fit
  integer :: i
  do i=1,9
    y(i)=1.5_dp+0.75_dp*(real(i-5,dp)/4.0_dp)**2
  end do
  theta0=[0.0_dp,0.0_dp]
  call nlrq_fit(y,theta0,0.5_dp,quad_model,fit)
  if (fit%info /= 0) error stop 'nlrq info'
  if (maxval(abs(fit%coefficients-[1.5_dp,0.75_dp])) > 1.0e-5_dp) error stop 'nlrq coef'
  print *, 'test_nonlinear: PASS'
contains
  subroutine quad_model(theta,fitted,jacobian)
    real(dp), intent(in) :: theta(:)
    real(dp), intent(out) :: fitted(:)
    real(dp), intent(out) :: jacobian(:,:)
    integer :: j
    real(dp) :: z
    do j=1,size(fitted)
      z=real(j-5,dp)/4.0_dp
      fitted(j)=theta(1)+theta(2)*z*z
      jacobian(j,1)=1.0_dp
      jacobian(j,2)=z*z
    end do
  end subroutine
end program
