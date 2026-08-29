program test_lbfgsb3
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use lbfgsb3_mod, only : dp, lbfgsb_control_t, lbfgsb_result_t, &
       lbfgsb_minimize, lbfgsb_minimize_fd, lbfgsb_success, &
       lbfgsb_iteration_limit, lbfgsb_input_error, &
       lbfgsb_evaluation_error, lbfgsb_user_stop
  implicit none


  type :: shift_data_t
    real(dp) :: target(3)
  end type shift_data_t

  integer :: failures
  failures = 0

  call test_bounds(failures)
  call test_rosenbrock(failures)
  call test_generalized_rosenbrock(failures)
  call test_finite_difference(failures)
  call test_chebyquad(failures)
  call test_user_data(failures)
  call test_limits(failures)
  call test_invalid_inputs(failures)
  call test_nonfinite(failures)
  call test_monitor_stop(failures)

  if (failures /= 0) then
    write(*,'(a,i0)') 'FAILED tests: ', failures
    error stop 1
  end if
  write(*,'(a)') 'All lbfgsb3 tests passed.'

contains

  subroutine test_bounds(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(4), lower(4), upper(4)
    type(lbfgsb_control_t) :: control
    type(lbfgsb_result_t) :: result
    integer :: i

    do i = 1, 4
      lower(i) = real((i-1)*3,dp)/4.0_dp
      upper(i) = real(i*5,dp)/4.0_dp
    end do
    x = 0.5_dp*(lower+upper)
    control%pgtol = 1.0e-10_dp
    control%reltol = 0.0_dp
    call lbfgsb_minimize(x, quadratic_fg, result, lower, upper, control)
    call check(maxval(abs(x - [0.0_dp,0.75_dp,1.5_dp,2.25_dp])) < 1.0e-10_dp, &
               'bound-constrained quadratic parameters', failures)
    call check(result%success, 'bound-constrained quadratic status', failures)
    call check(result%initial_point_projected .eqv. .false., &
               'feasible initial point not projected', failures)
  end subroutine test_bounds

  subroutine test_rosenbrock(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(2)
    type(lbfgsb_control_t) :: control
    type(lbfgsb_result_t) :: result

    x = [-1.2_dp, 1.0_dp]
    control%pgtol = 1.0e-9_dp
    control%reltol = 0.0_dp
    control%max_evaluations = 5000
    call lbfgsb_minimize(x, rosenbrock_fg, result, control=control)
    call check(maxval(abs(x - 1.0_dp)) < 2.0e-5_dp, &
               'Rosenbrock parameters', failures)
    call check(result%value < 1.0e-10_dp, 'Rosenbrock value', failures)
    call check(result%convergence == lbfgsb_success, &
               'Rosenbrock status', failures)
  end subroutine test_rosenbrock


  subroutine test_generalized_rosenbrock(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(100)
    type(lbfgsb_control_t) :: control
    type(lbfgsb_result_t) :: result

    x = 3.0_dp
    control%pgtol = 1.0e-8_dp
    control%reltol = 0.0_dp
    control%max_evaluations = 10000
    call lbfgsb_minimize(x, generalized_rosenbrock_fg, result, control=control)
    call check(abs(result%value-1.0_dp) < 1.0e-8_dp, &
               '100-parameter generalized Rosenbrock value', failures)
    call check(maxval(abs(x-1.0_dp)) < 2.0e-4_dp, &
               '100-parameter generalized Rosenbrock parameters', failures)
  end subroutine test_generalized_rosenbrock

  subroutine test_finite_difference(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(2), lower(1), upper(1)
    type(lbfgsb_control_t) :: control
    type(lbfgsb_result_t) :: result

    x = [-1.2_dp, 1.0_dp]
    lower = -2.0_dp
    upper = 2.0_dp
    control%pgtol = 1.0e-7_dp
    control%reltol = 0.0_dp
    control%max_evaluations = 20000
    call lbfgsb_minimize_fd(x, rosenbrock_f, result, lower, upper, control)
    call check(maxval(abs(x - 1.0_dp)) < 1.0e-4_dp, &
               'finite-difference Rosenbrock parameters', failures)
    call check(result%value < 1.0e-7_dp, &
               'finite-difference Rosenbrock value', failures)
    call check(result%function_evaluations > result%gradient_evaluations, &
               'finite-difference objective count', failures)
  end subroutine test_finite_difference

  subroutine test_chebyquad(failures)
    integer, intent(inout) :: failures
    integer, parameter :: ns(4) = [2,3,5,8]
    real(dp), parameter :: expected(8,4) = reshape([ &
      0.211324865405187_dp, 0.788675134594813_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.146446609406726_dp, 0.5_dp, 0.853553390593274_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.083751256499177_dp, 0.312729295223372_dp, 0.5_dp, 0.687270704776628_dp, 0.916248743500823_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.043_dp, 0.193_dp, 0.266_dp, 0.5_dp, 0.5_dp, 0.734_dp, 0.807_dp, 0.957_dp ], [8,4])
    real(dp), allocatable :: x(:), lower(:), upper(:)
    type(lbfgsb_control_t) :: control
    type(lbfgsb_result_t) :: result
    integer :: k, n, i
    real(dp) :: tol

    control%pgtol = 1.0e-10_dp
    control%reltol = 0.0_dp
    control%max_evaluations = 10000
    do k = 1, size(ns)
      n = ns(k)
      allocate(x(n), lower(n), upper(n))
      do i = 1, n
        x(i) = real(i,dp)/real(n+1,dp)
      end do
      lower = -10.0_dp
      upper = 10.0_dp
      call lbfgsb_minimize(x, chebyquad_fg, result, lower, upper, control)
      tol = merge(8.0e-3_dp, 8.0e-5_dp, n == 8)
      call check(maxval(abs(x - expected(1:n,k))) < tol, &
                 'Chebyquad parameters n='//integer_string(n), failures)
      if (n < 8) then
        call check(result%value < 1.0e-10_dp, &
                   'Chebyquad value n='//integer_string(n), failures)
      else
        call check(abs(result%value - 3.5169e-3_dp) < 8.0e-4_dp, &
                   'Chebyquad value n=8', failures)
      end if
      deallocate(x, lower, upper)
    end do
  end subroutine test_chebyquad

  subroutine test_user_data(failures)
    integer, intent(inout) :: failures
    type(shift_data_t) :: data
    real(dp) :: x(3)
    type(lbfgsb_control_t) :: control
    type(lbfgsb_result_t) :: result

    data%target = [2.0_dp, -1.0_dp, 0.5_dp]
    x = 0.0_dp
    control%pgtol = 1.0e-12_dp
    control%reltol = 0.0_dp
    call lbfgsb_minimize(x, shifted_fg, result, control=control, user_data=data)
    call check(maxval(abs(x-data%target)) < 1.0e-9_dp, &
               'polymorphic user data', failures)
  end subroutine test_user_data

  subroutine test_limits(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(2)
    type(lbfgsb_control_t) :: control
    type(lbfgsb_result_t) :: result

    x = [-1.2_dp, 1.0_dp]
    control%max_evaluations = 1
    control%reltol = 0.0_dp
    call lbfgsb_minimize(x, rosenbrock_fg, result, control=control)
    call check(result%convergence == lbfgsb_iteration_limit, &
               'evaluation limit status', failures)
  end subroutine test_limits

  subroutine test_invalid_inputs(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(2), lower(2), upper(2)
    type(lbfgsb_result_t) :: result

    x = 0.0_dp
    lower = [1.0_dp, -1.0_dp]
    upper = [0.0_dp, 1.0_dp]
    call lbfgsb_minimize(x, quadratic_fg, result, lower, upper)
    call check(result%convergence == lbfgsb_input_error, &
               'invalid bounds status', failures)
    call check(result%offending_index == 1, &
               'invalid bounds index', failures)
  end subroutine test_invalid_inputs

  subroutine test_nonfinite(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(1)
    type(lbfgsb_result_t) :: result

    x = 0.0_dp
    call lbfgsb_minimize(x, nan_fg, result)
    call check(result%convergence == lbfgsb_evaluation_error, &
               'non-finite callback status', failures)
  end subroutine test_nonfinite

  subroutine test_monitor_stop(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(2)
    type(lbfgsb_control_t) :: control
    type(lbfgsb_result_t) :: result

    x = [-1.2_dp, 1.0_dp]
    control%reltol = 0.0_dp
    call lbfgsb_minimize(x, rosenbrock_fg, result, control=control, monitor=stop_after_two)
    call check(result%convergence == lbfgsb_user_stop, &
               'monitor cancellation status', failures)
    call check(result%iterations == 2, &
               'monitor cancellation iteration', failures)
  end subroutine test_monitor_stop

  subroutine quadratic_fg(x, f, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    class(*), intent(inout), optional :: user_data
    f = dot_product(x,x)
    g = 2.0_dp*x
  end subroutine quadratic_fg

  subroutine rosenbrock_fg(x, f, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    class(*), intent(inout), optional :: user_data
    f = 100.0_dp*(x(2)-x(1)*x(1))**2 + (1.0_dp-x(1))**2
    g(1) = -400.0_dp*x(1)*(x(2)-x(1)*x(1)) - 2.0_dp*(1.0_dp-x(1))
    g(2) = 200.0_dp*(x(2)-x(1)*x(1))
  end subroutine rosenbrock_fg


  subroutine generalized_rosenbrock_fg(x, f, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    class(*), intent(inout), optional :: user_data
    integer :: n, i
    real(dp), parameter :: scale = 10.0_dp
    real(dp) :: z

    n = size(x)
    f = 1.0_dp
    g = 0.0_dp
    do i = 1, n-1
      z = x(i)*x(i)-x(i+1)
      f = f+scale*z*z+(x(i+1)-1.0_dp)**2
      g(i) = g(i)+4.0_dp*scale*x(i)*z
      g(i+1) = g(i+1)-2.0_dp*scale*z+2.0_dp*(x(i+1)-1.0_dp)
    end do
  end subroutine generalized_rosenbrock_fg

  function rosenbrock_f(x, user_data) result(f)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: f
    f = 100.0_dp*(x(2)-x(1)*x(1))**2 + (1.0_dp-x(1))**2
  end function rosenbrock_f

  subroutine chebyquad_fg(x, f, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    class(*), intent(inout), optional :: user_data
    real(dp), allocatable :: residual(:), jacobian(:,:)
    integer :: n, i, k, j
    real(dp) :: z2, z4, z5, z6, z7, z8, rr

    n = size(x)
    allocate(residual(n), jacobian(n,n))
    residual = 0.0_dp
    jacobian = 0.0_dp
    do i = 1, n
      rr = 0.0_dp
      do k = 1, n
        z7 = 1.0_dp
        z2 = 2.0_dp*x(k)-1.0_dp
        z8 = z2
        do j = 1, i-1
          z6 = z7
          z7 = z8
          z8 = 2.0_dp*z2*z7-z6
        end do
        rr = rr+z8
      end do
      rr = rr/real(n,dp)
      if (mod(i,2) == 0) rr = rr+1.0_dp/real(i*i-1,dp)
      residual(i) = rr

      do k = 1, n
        z5 = 0.0_dp
        jacobian(i,k) = 2.0_dp
        z8 = 2.0_dp*x(k)-1.0_dp
        z2 = z8
        z7 = 1.0_dp
        do j = 1, i-1
          z4 = z5
          z5 = jacobian(i,k)
          jacobian(i,k) = 4.0_dp*z8+2.0_dp*z2*z5-z4
          z6 = z7
          z7 = z8
          z8 = 2.0_dp*z2*z7-z6
        end do
        jacobian(i,k) = jacobian(i,k)/real(n,dp)
      end do
    end do
    f = dot_product(residual,residual)
    g = 2.0_dp*matmul(transpose(jacobian), residual)
  end subroutine chebyquad_fg

  subroutine shifted_fg(x, f, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    class(*), intent(inout), optional :: user_data
    select type (data => user_data)
    type is (shift_data_t)
      f = sum((x-data%target)**2)
      g = 2.0_dp*(x-data%target)
    class default
      f = huge(1.0_dp)
      g = 0.0_dp
    end select
  end subroutine shifted_fg

  subroutine nan_fg(x, f, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    class(*), intent(inout), optional :: user_data
    f = ieee_value(0.0_dp, ieee_quiet_nan)
    g = 0.0_dp
  end subroutine nan_fg

  subroutine stop_after_two(x, f, g, iteration, evaluations, stop, user_data)
    real(dp), intent(in) :: x(:), f, g(:)
    integer, intent(in) :: iteration, evaluations
    logical, intent(out) :: stop
    class(*), intent(inout), optional :: user_data
    stop = iteration >= 2
  end subroutine stop_after_two

  subroutine check(condition, name, failures)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: name
    integer, intent(inout) :: failures
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(name)
      failures = failures+1
    end if
  end subroutine check

  function integer_string(i) result(text)
    integer, intent(in) :: i
    character(len=:), allocatable :: text
    character(len=32) :: buffer
    write(buffer,'(i0)') i
    text = trim(buffer)
  end function integer_string

end program test_lbfgsb3
