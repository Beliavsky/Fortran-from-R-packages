program nonlinear_least_squares
  use optimflex
  implicit none
  type(optim_result) :: res
  real(dp) :: x0(2)
  x0=[0.0_dp,0.0_dp]
  call gauss_newton(x0,obj,res,residual=residuals)
  print '(a,2f14.8,a,es12.4)', 'Gauss-Newton par = ',res%par,'  f = ',res%objective
contains
  real(dp) function obj(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: r(:)
    r=residuals(x); f=sum(r*r)
  end function obj
  function residuals(x) result(r)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: r(:)
    allocate(r(2)); r=[x(1)-2.0_dp,x(2)+1.0_dp]
  end function residuals
end program nonlinear_least_squares
