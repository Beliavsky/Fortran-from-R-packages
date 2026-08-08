program test_projection
  use bb, only: dp, spg_linear, spg_control, spg_result
  implicit none
  real(dp) :: p0(2), a(4,2), b(4)
  type(spg_control) :: ctrl
  type(spg_result) :: fit

  p0 = [-0.3599199_dp, -1.2219309_dp]
  a(1,:) = [ 1.0_dp, -1.0_dp]
  a(2,:) = [ 1.0_dp,  1.0_dp]
  a(3,:) = [-1.0_dp,  1.0_dp]
  a(4,:) = [-1.0_dp, -1.0_dp]
  b = -1.0_dp
  ctrl = spg_control(method=3, m=10, maxit=3000, gtol=1.0e-7_dp, ftol=1.0e-13_dp)

  fit = spg_linear(p0, objective, a, b, 0, ctrl, gradient)
  if (.not. fit%succeeded()) error stop 'linear projection inequality solve failed'
  if (maxval(abs(fit%par - [1.0_dp, 0.0_dp])) > 2.0e-5_dp) &
    error stop 'linear projection inequality result wrong'

  fit = spg_linear(p0, objective, a, b, 1, ctrl, gradient)
  if (.not. fit%succeeded()) error stop 'linear projection equality solve failed'
  if (maxval(abs(fit%par - [0.0_dp, 1.0_dp])) > 2.0e-5_dp) &
    error stop 'linear projection equality result wrong'

  a = 0.0_dp
  a(1,:) = [ 1.0_dp,  0.0_dp]
  a(2,:) = [ 0.0_dp,  1.0_dp]
  a(3,:) = [-1.0_dp,  0.0_dp]
  a(4,:) = [ 0.0_dp, -1.0_dp]
  b = [0.0_dp, 0.0_dp, -0.5_dp, -0.5_dp]
  fit = spg_linear(p0, objective, a, b, 0, ctrl, gradient)
  if (.not. fit%succeeded()) error stop 'linear projection box-equivalent solve failed'
  if (abs(fit%par(1)-0.5_dp) > 2.0e-6_dp .or. abs(fit%par(2)-0.125_dp) > 2.0e-2_dp) &
    error stop 'linear projection box-equivalent result wrong'

  print *, 'test_projection: PASS'
contains
  function objective(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = (x(1)-1.5_dp)**2 + (x(2)-0.125_dp)**4
  end function objective

  subroutine gradient(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    g(1) = 2.0_dp*(x(1)-1.5_dp)
    g(2) = 4.0_dp*(x(2)-0.125_dp)**3
  end subroutine gradient
end program test_projection
