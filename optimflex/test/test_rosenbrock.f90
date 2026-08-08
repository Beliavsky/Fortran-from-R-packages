program test_rosenbrock
  use optimflex
  implicit none
  type(optim_result) :: r
  type(optim_control) :: c
  real(dp) :: x0(2), lo(2), hi(2)
  x0 = [-1.2_dp, 1.0_dp]
  lo = -5.0_dp
  hi = 5.0_dp

  c = bfgs_default_control(); c%max_iter = 5000
  call bfgs(x0,rosen,r,control=c); call check(r,'bfgs')
  c = lbfgsb_default_control(); c%max_iter = 5000
  call l_bfgs_b(x0,rosen,r,lo,hi,control=c); call check(r,'l_bfgs_b')
  c = newton_default_control(); c%max_iter = 5000; c%diff_method = diff_central
  call newton_raphson(x0,rosen,r,control=c); call check(r,'newton_raphson')
  c = modified_newton_default_control(); c%max_iter = 5000; c%diff_method = diff_central
  call modified_newton(x0,rosen,r,control=c); call check(r,'modified_newton')
  c = dogleg_default_control(); c%max_iter = 5000; c%initial_delta = 2.0_dp
  call dogleg(x0,rosen,r,lower=lo,upper=hi,control=c); call check(r,'dogleg')
  c = double_dogleg_default_control(); c%max_iter = 5000; c%initial_delta = 2.0_dp
  call double_dogleg(x0,rosen,r,lower=lo,upper=hi,control=c); call check(r,'double_dogleg')
  c = lm_default_control(); c%max_iter = 5000
  call levenberg_marquardt(x0,rosen,r,lower=lo,upper=hi,control=c); call check(r,'levenberg_marquardt')
  print *, 'test_rosenbrock: PASS'
contains
  real(dp) function rosen(x) result(f)
    real(dp), intent(in) :: x(:)
    f = 100.0_dp*(x(2)-x(1)*x(1))**2 + (1.0_dp-x(1))**2
  end function rosen
  subroutine check(res,name)
    type(optim_result), intent(in) :: res
    character(len=*), intent(in) :: name
    if (.not. res%converged) then
      print *, trim(name), trim(res%status), res%par, res%objective
      error stop 'Rosenbrock convergence failed'
    end if
    if (maxval(abs(res%par-[1.0_dp,1.0_dp])) > 2.0e-2_dp) error stop 'Rosenbrock parameter mismatch'
  end subroutine check
end program test_rosenbrock
