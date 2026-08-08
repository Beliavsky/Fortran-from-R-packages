program test_local
  use quantreg, only : dp, lprq_result, lprq
  implicit none
  real(dp) :: x(11), y(11)
  type(lprq_result) :: fit
  integer :: i
  do i=1,11
    x(i)=real(i-6,dp)/5.0_dp
  end do
  y=1.0_dp+2.0_dp*x
  call lprq(x,y,0.6_dp,0.5_dp,7,fit)
  if (fit%info /= 0) error stop 'lprq info'
  if (maxval(abs(fit%fitted-(1.0_dp+2.0_dp*fit%x))) > 1.0e-5_dp) error stop 'lprq fit'
  if (maxval(abs(fit%derivative-2.0_dp)) > 1.0e-5_dp) error stop 'lprq slope'
  print *, 'test_local: PASS'
end program
