program test_gradient
  use calibrar, only : dp, numerical_gradient, quadratic_objective
  implicit none
  real(dp) :: x(3), g(3)
  character(len=10), parameter :: methods(4)=[character(len=10)::"forward","backward","central","richardson"]
  integer :: i
  x=[1.0_dp,-2.0_dp,0.5_dp]
  do i=1,4
    call numerical_gradient(quadratic_objective,x,g,trim(methods(i)))
    if(maxval(abs(g-2.0_dp*x))>2.0e-5_dp) error stop "gradient test failed"
  end do
  print *, "PASS test_gradient"
end program test_gradient
