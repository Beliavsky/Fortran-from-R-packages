program test_drivers
  use bb, only: dp, bboptim, bbsolve, multistart_optimize, spg_result, sane_result, &
    multistart_result, bboptim_control, bbsolve_control
  implicit none
  type(spg_result) :: opt
  type(sane_result) :: sol
  type(multistart_result) :: ms
  type(bboptim_control) :: oc
  type(bbsolve_control) :: sc
  real(dp) :: starts(2,3)

  oc = bboptim_control()
  opt = bboptim([5.0_dp, -5.0_dp], quad, oc, quad_grad)
  if (.not. opt%succeeded()) error stop 'BBoptim failed'
  if (maxval(abs(opt%par - [2.0_dp, -1.0_dp])) > 2.0e-5_dp) error stop 'BBoptim wrong'

  sc = bbsolve_control(try_nm=.false.)
  sol = bbsolve([3.0_dp, 3.0_dp], system_fn, sc)
  if (.not. sol%succeeded()) error stop 'BBsolve failed'
  if (maxval(abs(sol%par - [1.0_dp, -2.0_dp])) > 2.0e-5_dp) error stop 'BBsolve wrong'

  starts(:,1) = [-4.0_dp, 4.0_dp]
  starts(:,2) = [0.0_dp, 0.0_dp]
  starts(:,3) = [6.0_dp, -8.0_dp]
  ms = multistart_optimize(starts, quad, oc, quad_grad)
  if (.not. all(ms%converged)) error stop 'multistart optimize failed'
  if (maxval(abs(ms%par - spread([2.0_dp,-1.0_dp],2,3))) > 3.0e-5_dp) &
    error stop 'multistart optimize wrong'

  print *, 'test_drivers: PASS'
contains
  function quad(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = (x(1)-2.0_dp)**2 + 3.0_dp*(x(2)+1.0_dp)**2
  end function quad
  subroutine quad_grad(x,g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    g = [2.0_dp*(x(1)-2.0_dp), 6.0_dp*(x(2)+1.0_dp)]
  end subroutine quad_grad
  subroutine system_fn(x,f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f = [x(1)-1.0_dp, x(2)+2.0_dp]
  end subroutine system_fn
end program test_drivers
