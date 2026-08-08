program rosenbrock_example
  use mize_mod, only : dp, mize_control_t, mize_result_t, mize_minimize
  implicit none

  type(mize_control_t) :: control
  type(mize_result_t) :: result
  real(dp) :: x(2)

  x = [-1.2_dp, 1.0_dp]
  control = mize_control_t()
  control%method = 'L-BFGS'
  control%memory = 7
  control%max_iterations = 500
  control%grad_tol = 1.0e-8_dp
  control%ginf_tol = 1.0e-8_dp
  control%store_progress = .true.

  call mize_minimize(x, rosenbrock_fg, result, control)

  write(*, '(a, 2(1x, es16.8))') 'parameters:', x
  write(*, '(a, 1x, es16.8)') 'objective:', result%value
  write(*, '(a, 1x, i0)') 'iterations:', result%iterations
  write(*, '(a, 1x, i0)') 'evaluations:', result%function_evaluations
  write(*, '(a, 1x, a)') 'status:', result%message

contains

  subroutine rosenbrock_fg(x, f, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f
    real(dp), intent(out) :: g(:)
    class(*), intent(inout), optional :: user_data

    f = 100.0_dp * (x(2) - x(1) * x(1)) ** 2 + (1.0_dp - x(1)) ** 2
    g(1) = -400.0_dp * x(1) * (x(2) - x(1) * x(1)) - 2.0_dp * (1.0_dp - x(1))
    g(2) = 200.0_dp * (x(2) - x(1) * x(1))
  end subroutine rosenbrock_fg

end program rosenbrock_example
