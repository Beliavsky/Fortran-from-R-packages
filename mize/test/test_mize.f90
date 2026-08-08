program test_mize
  use mize_mod, only : dp, mize_control_t, mize_result_t, mize_state_t, &
       gradient_check_t, mize_minimize, mize_init, mize_step, &
       mize_state_result, check_mize_gradient, mize_success, &
       mize_user_stop
  implicit none

  type :: scale_data_t
    real(dp) :: scale = 2.0_dp
  end type scale_data_t

  integer :: failures

  failures = 0
  call test_quasi_newton(failures)
  call test_cg_updates(failures)
  call test_newton_methods(failures)
  call test_first_order_methods(failures)
  call test_line_search_aliases(failures)
  call test_stateful_api(failures)
  call test_gradient_diagnostic(failures)
  call test_user_data_and_monitor(failures)

  if (failures /= 0) then
    write(*, '(a, i0)') 'FAILED tests: ', failures
    error stop 1
  end if
  write(*, '(a)') 'All mize tests passed.'

contains

  subroutine test_quasi_newton(failures)
    integer, intent(inout) :: failures
    character(len=8), parameter :: methods(3) = [character(len=8) :: &
         'BFGS', 'SR1', 'L-BFGS']
    type(mize_control_t) :: control
    type(mize_result_t) :: result
    real(dp) :: x(2)
    integer :: i

    do i = 1, size(methods)
      x = [-1.2_dp, 1.0_dp]
      control = mize_control_t()
      control%method = methods(i)
      control%max_iterations = 500
      control%grad_tol = 1.0e-7_dp
      control%ginf_tol = 1.0e-7_dp
      control%abs_tol = 1.0e-14_dp
      control%rel_tol = 1.0e-14_dp
      control%step_tol = 1.0e-14_dp
      call mize_minimize(x, rosenbrock_fg, result, control)
      call assert_true(result%status == mize_success, trim(methods(i))//' status', failures)
      call assert_close(x(1), 1.0_dp, 2.0e-4_dp, trim(methods(i))//' x1', failures)
      call assert_close(x(2), 1.0_dp, 2.0e-4_dp, trim(methods(i))//' x2', failures)
      call assert_true(result%value < 1.0e-8_dp, trim(methods(i))//' value', failures)
    end do
  end subroutine test_quasi_newton

  subroutine test_cg_updates(failures)
    integer, intent(inout) :: failures
    character(len=4), parameter :: updates(6) = [character(len=4) :: &
         'FR', 'PR+', 'HS+', 'DY', 'HZ+', 'PRFR']
    type(mize_control_t) :: control
    type(mize_result_t) :: result
    real(dp) :: x(4)
    integer :: i

    do i = 1, size(updates)
      x = [4.0_dp, -3.0_dp, 2.0_dp, -1.0_dp]
      control = mize_control_t()
      control%method = 'CG'
      control%cg_update = updates(i)
      control%c2 = 0.1_dp
      control%max_iterations = 300
      control%grad_tol = 1.0e-8_dp
      control%ginf_tol = 1.0e-8_dp
      control%abs_tol = 1.0e-14_dp
      control%rel_tol = 1.0e-14_dp
      call mize_minimize(x, quadratic_fg, result, control)
      call assert_true(result%status == mize_success, 'CG '//trim(updates(i))//' status', failures)
      call assert_true(maxval(abs(x)) < 1.0e-5_dp, 'CG '//trim(updates(i))//' solution', failures)
    end do
  end subroutine test_cg_updates

  subroutine test_newton_methods(failures)
    integer, intent(inout) :: failures
    character(len=6), parameter :: methods(3) = [character(len=6) :: &
         'Newton', 'pHess', 'TN']
    type(mize_control_t) :: control
    type(mize_result_t) :: result
    real(dp) :: x(4)
    integer :: i

    do i = 1, size(methods)
      x = [3.0_dp, -2.0_dp, 1.0_dp, 4.0_dp]
      control = mize_control_t()
      control%method = methods(i)
      control%max_iterations = 100
      control%grad_tol = 1.0e-10_dp
      control%ginf_tol = 1.0e-10_dp
      control%abs_tol = 1.0e-14_dp
      control%rel_tol = 1.0e-14_dp
      if (trim(methods(i)) == 'TN') then
        call mize_minimize(x, quadratic_fg, result, control, hvp=quadratic_hvp)
      else
        call mize_minimize(x, quadratic_fg, result, control, hessian=quadratic_hessian)
      end if
      call assert_true(result%status == mize_success, trim(methods(i))//' status', failures)
      call assert_true(maxval(abs(x)) < 1.0e-8_dp, trim(methods(i))//' solution', failures)
    end do
  end subroutine test_newton_methods

  subroutine test_first_order_methods(failures)
    integer, intent(inout) :: failures
    character(len=8), parameter :: methods(4) = [character(len=8) :: &
         'SD', 'MOM', 'NAG', 'DBD']
    type(mize_control_t) :: control
    type(mize_result_t) :: result
    real(dp) :: x(3)
    integer :: i

    do i = 1, size(methods)
      x = [2.0_dp, -1.5_dp, 1.0_dp]
      control = mize_control_t()
      control%method = methods(i)
      control%line_search = 'constant'
      control%step0 = 0.05_dp
      control%mom_init = 0.5_dp
      control%max_iterations = 1000
      control%grad_tol = 1.0e-6_dp
      control%ginf_tol = 1.0e-6_dp
      control%abs_tol = -1.0_dp
      control%rel_tol = -1.0_dp
      control%step_tol = -1.0_dp
      if (trim(methods(i)) == 'DBD') then
        control%step0 = 0.01_dp
        control%step_up = 1.02_dp
        control%step_down = 0.5_dp
      end if
      call mize_minimize(x, identity_quadratic_fg, result, control)
      call assert_true(result%status == mize_success, trim(methods(i))//' status', failures)
      call assert_true(maxval(abs(x)) < 2.0e-4_dp, trim(methods(i))//' solution', failures)
    end do
  end subroutine test_first_order_methods

  subroutine test_line_search_aliases(failures)
    integer, intent(inout) :: failures
    character(len=14), parameter :: searches(7) = [character(len=14) :: &
         'More-Thuente', 'Rasmussen', 'Schmidt', 'Hager-Zhang', &
         'backtracking', 'bold driver', 'constant']
    type(mize_control_t) :: control
    type(mize_result_t) :: result
    real(dp) :: x(2)
    integer :: i

    do i = 1, size(searches)
      x = [2.0_dp, -1.0_dp]
      control = mize_control_t()
      control%method = 'BFGS'
      control%line_search = searches(i)
      control%step0 = merge(0.1_dp, 1.0_dp, trim(searches(i)) == 'constant')
      control%max_iterations = 300
      control%grad_tol = 1.0e-7_dp
      control%ginf_tol = 1.0e-7_dp
      control%abs_tol = 1.0e-14_dp
      control%rel_tol = 1.0e-14_dp
      call mize_minimize(x, identity_quadratic_fg, result, control)
      call assert_true(result%status == mize_success, trim(searches(i))//' status', failures)
      call assert_true(maxval(abs(x)) < 1.0e-4_dp, trim(searches(i))//' solution', failures)
    end do
  end subroutine test_line_search_aliases

  subroutine test_stateful_api(failures)
    integer, intent(inout) :: failures
    type(mize_control_t) :: control
    type(mize_state_t) :: state
    type(mize_result_t) :: result
    real(dp) :: x(2)

    x = [-1.2_dp, 1.0_dp]
    control = mize_control_t()
    control%method = 'L-BFGS'
    control%max_iterations = 500
    control%grad_tol = 1.0e-7_dp
    control%ginf_tol = 1.0e-7_dp
    control%store_progress = .true.
    call mize_init(state, x, rosenbrock_fg, control)
    do while (.not. state%terminated)
      call mize_step(state, rosenbrock_fg)
    end do
    call mize_state_result(state, result)
    call assert_true(result%status == mize_success, 'stateful status', failures)
    call assert_true(result%iterations > 0, 'stateful iterations', failures)
    call assert_true(allocated(result%progress_value), 'stateful progress', failures)
    call assert_true(result%value < 1.0e-8_dp, 'stateful value', failures)
  end subroutine test_stateful_api

  subroutine test_gradient_diagnostic(failures)
    integer, intent(inout) :: failures
    type(gradient_check_t) :: check
    real(dp) :: x(2)

    x = [-0.8_dp, 1.2_dp]
    call check_mize_gradient(x, rosenbrock_fg, check, tolerance=1.0e-5_dp)
    call assert_true(check%passed, 'gradient check passed', failures)
    call assert_true(check%maximum_relative_error < 1.0e-5_dp, &
                     'gradient check error', failures)
  end subroutine test_gradient_diagnostic

  subroutine test_user_data_and_monitor(failures)
    integer, intent(inout) :: failures
    type(scale_data_t) :: data
    type(mize_control_t) :: control
    type(mize_result_t) :: result
    real(dp) :: x(2)

    x = [2.0_dp, -1.0_dp]
    control = mize_control_t()
    control%method = 'SD'
    control%line_search = 'constant'
    control%step0 = 0.1_dp
    control%max_iterations = 100
    control%grad_tol = -1.0_dp
    control%ginf_tol = -1.0_dp
    control%abs_tol = -1.0_dp
    control%rel_tol = -1.0_dp
    control%step_tol = -1.0_dp
    call mize_minimize(x, scaled_quadratic_fg, result, control, &
                       monitor=stop_monitor, user_data=data)
    call assert_true(result%status == mize_user_stop, 'monitor status', failures)
    call assert_true(result%iterations == 3, 'monitor iteration', failures)
  end subroutine test_user_data_and_monitor

  subroutine rosenbrock_fg(x, f, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    class(*), intent(inout), optional :: user_data
    f = 100.0_dp * (x(2) - x(1) * x(1)) ** 2 + (1.0_dp - x(1)) ** 2
    g(1) = -400.0_dp * x(1) * (x(2) - x(1) * x(1)) - 2.0_dp * (1.0_dp - x(1))
    g(2) = 200.0_dp * (x(2) - x(1) * x(1))
  end subroutine rosenbrock_fg

  subroutine quadratic_fg(x, f, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    class(*), intent(inout), optional :: user_data
    integer :: i
    f = 0.0_dp
    do i = 1, size(x)
      f = f + 0.5_dp * real(i, dp) * x(i) * x(i)
      g(i) = real(i, dp) * x(i)
    end do
  end subroutine quadratic_fg

  subroutine identity_quadratic_fg(x, f, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    class(*), intent(inout), optional :: user_data
    f = 0.5_dp * dot_product(x, x)
    g = x
  end subroutine identity_quadratic_fg

  subroutine quadratic_hessian(x, h, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: h(:, :)
    class(*), intent(inout), optional :: user_data
    integer :: i
    h = 0.0_dp
    do i = 1, size(x)
      h(i, i) = real(i, dp)
    end do
  end subroutine quadratic_hessian

  subroutine quadratic_hvp(x, v, hv, user_data)
    real(dp), intent(in) :: x(:), v(:)
    real(dp), intent(out) :: hv(:)
    class(*), intent(inout), optional :: user_data
    integer :: i
    do i = 1, size(x)
      hv(i) = real(i, dp) * v(i)
    end do
  end subroutine quadratic_hvp

  subroutine scaled_quadratic_fg(x, f, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: scale
    scale = 1.0_dp
    if (present(user_data)) then
      select type (user_data)
      type is (scale_data_t)
        scale = user_data%scale
      end select
    end if
    f = 0.5_dp * scale * dot_product(x, x)
    g = scale * x
  end subroutine scaled_quadratic_fg

  subroutine stop_monitor(x, f, g, iteration, function_evaluations, &
                          gradient_evaluations, stop, user_data)
    real(dp), intent(in) :: x(:), f, g(:)
    integer, intent(in) :: iteration, function_evaluations, gradient_evaluations
    logical, intent(out) :: stop
    class(*), intent(inout), optional :: user_data
    stop = iteration >= 3
  end subroutine stop_monitor

  subroutine assert_true(condition, label, failures)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failures
    if (.not. condition) then
      write(*, '(a)') 'FAIL: '//trim(label)
      failures = failures + 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tolerance, label, failures)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failures
    call assert_true(abs(actual - expected) <= tolerance, label, failures)
  end subroutine assert_close

end program test_mize
