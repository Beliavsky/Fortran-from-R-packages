program spg_rosenbrock
  use bb, only: dp, spg, spg_control, spg_result
  implicit none

  type(spg_control) :: control
  type(spg_result) :: fit

  control = spg_control(method=3, m=10, maxit=5000, gtol=1.0e-7_dp)
  fit = spg([-1.2_dp, 1.0_dp], objective, control, gradient)

  print '(a,2f14.8)', 'par = ', fit%par
  print '(a,es14.6)', 'value = ', fit%value
  print '(a,i0)', 'convergence = ', fit%convergence
  print '(a,a)', 'message = ', fit%message

contains

  function objective(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = 100.0_dp * (x(2) - x(1)**2)**2 + (1.0_dp - x(1))**2
  end function objective

  subroutine gradient(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    g(1) = -400.0_dp*x(1)*(x(2)-x(1)**2) - 2.0_dp*(1.0_dp-x(1))
    g(2) = 200.0_dp*(x(2)-x(1)**2)
  end subroutine gradient

end program spg_rosenbrock
