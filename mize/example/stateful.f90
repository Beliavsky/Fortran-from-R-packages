program stateful_example
  use mize_mod, only : dp, mize_control_t, mize_state_t, mize_result_t, &
       mize_init, mize_step, mize_state_result
  implicit none

  type(mize_control_t) :: control
  type(mize_state_t) :: optimizer
  type(mize_result_t) :: result
  real(dp) :: x(3)

  x = [4.0_dp, -3.0_dp, 2.0_dp]
  control = mize_control_t()
  control%method = 'CG'
  control%cg_update = 'HZ+'
  control%c2 = 0.1_dp
  control%max_iterations = 100

  call mize_init(optimizer, x, quadratic_fg, control)
  do while (.not. optimizer%terminated)
    call mize_step(optimizer, quadratic_fg)
  end do
  call mize_state_result(optimizer, result)

  write(*, '(a, 3(1x, es14.6))') 'parameters:', result%par
  write(*, '(a, 1x, es14.6)') 'objective:', result%value
  write(*, '(a, 1x, a)') 'status:', result%message

contains

  subroutine quadratic_fg(x, f, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f
    real(dp), intent(out) :: g(:)
    class(*), intent(inout), optional :: user_data
    integer :: i

    f = 0.0_dp
    do i = 1, size(x)
      f = f + 0.5_dp * real(i, dp) * x(i) * x(i)
      g(i) = real(i, dp) * x(i)
    end do
  end subroutine quadratic_fg

end program stateful_example
