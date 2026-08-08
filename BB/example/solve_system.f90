program solve_system
  use bb, only: dp, bbsolve, sane_result
  implicit none

  type(sane_result) :: fit

  fit = bbsolve([4.0_dp, -5.0_dp], equations)

  print '(a,2f14.8)', 'root = ', fit%par
  print '(a,es14.6)', 'residual = ', fit%residual
  print '(a,i0)', 'convergence = ', fit%convergence
  print '(a,a)', 'message = ', fit%message

contains

  subroutine equations(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(1) = x(1) - 1.0_dp
    f(2) = x(2) + 2.0_dp
  end subroutine equations

end program solve_system
