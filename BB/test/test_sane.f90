program test_sane
  use bb, only: dp, sane, sane_control, sane_result
  implicit none
  type(sane_control) :: ctrl
  type(sane_result) :: fit

  ctrl = sane_control(method=2, m=10, maxit=1000, tol=1.0e-9_dp, trace=.false.)
  fit = sane([4.0_dp, -5.0_dp], system_fn, ctrl)
  if (.not. fit%succeeded()) error stop 'SANE did not converge'
  if (maxval(abs(fit%par - [1.0_dp, -2.0_dp])) > 2.0e-6_dp) error stop 'SANE root wrong'
  if (fit%residual > 1.0e-8_dp) error stop 'SANE residual too large'
  print *, 'test_sane: PASS'
contains
  subroutine system_fn(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(1) = x(1) - 1.0_dp
    f(2) = x(2) + 2.0_dp
  end subroutine system_fn
end program test_sane
