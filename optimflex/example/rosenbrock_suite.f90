program rosenbrock_suite
  use optimflex
  implicit none
  type(optim_result) :: res
  real(dp) :: x0(2)
  x0=[-1.2_dp,1.0_dp]
  call bfgs(x0,rosen,res)
  print '(a,2f14.8,a,es12.4)', 'BFGS par = ',res%par,'  f = ',res%objective
contains
  real(dp) function rosen(x) result(f)
    real(dp), intent(in) :: x(:)
    f=100.0_dp*(x(2)-x(1)*x(1))**2+(1.0_dp-x(1))**2
  end function rosen
end program rosenbrock_suite
