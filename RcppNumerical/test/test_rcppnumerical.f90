program test_rcppnumerical
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
  use rcppnumerical
  implicit none

  integer :: failures
  type shift_data_t
    real(dp) :: target = 0.0_dp
  end type shift_data_t

  failures = 0
  call test_gauss_kronrod(failures)
  call test_infinite_1d(failures)
  call test_cuhre(failures)
  call test_optimization(failures)
  call test_fast_lr(failures)

  if (failures /= 0) then
    write(*,'(a,i0)') 'FAILED tests: ', failures
    error stop 1
  end if
  write(*,'(a)') 'All RcppNumerical Fortran tests passed.'

contains

  subroutine test_gauss_kronrod(failures)
    integer, intent(inout) :: failures
    type(integration_result_t) :: result
    integer :: rule

    do rule = gk15, gk201
      call integrate_1d(power_six, 0.0_dp, 1.0_dp, result, &
                        eps_abs=1.0e-12_dp, eps_rel=1.0e-12_dp, rule=rule)
      call check_close('Gauss-Kronrod rule', result%value, 1.0_dp/7.0_dp, &
                       2.0e-12_dp, failures)
      call check_true('Gauss-Kronrod status', result%error_code == 0, failures)
    end do

    call integrate_1d(beta_pdf, 0.3_dp, 0.8_dp, result)
    call check_close('beta probability', result%value, &
                     0.2528108217749998_dp, 2.0e-12_dp, failures)
    call integrate_1d(power_six, 1.0_dp, 0.0_dp, result)
    call check_close('reversed interval', result%value, -1.0_dp/7.0_dp, &
                     2.0e-12_dp, failures)
    call integrate_1d(power_six, 0.0_dp, 1.0_dp, result, rule=99)
    call check_true('invalid rule', result%error_code == 6, failures)
  end subroutine test_gauss_kronrod

  subroutine test_infinite_1d(failures)
    integer, intent(inout) :: failures
    type(integration_result_t) :: result
    real(dp) :: inf

    inf = ieee_value(1.0_dp, ieee_positive_inf)
    call integrate_1d(gaussian_kernel, -inf, inf, result, &
                      eps_abs=1.0e-11_dp, eps_rel=1.0e-11_dp)
    call check_close('doubly infinite integral', result%value, &
                     sqrt(acos(-1.0_dp)), 2.0e-10_dp, failures)
    call integrate_1d(exponential_tail, 0.0_dp, inf, result, &
                      eps_abs=1.0e-11_dp, eps_rel=1.0e-11_dp)
    call check_close('semi-infinite integral', result%value, 1.0_dp, &
                     2.0e-10_dp, failures)
  end subroutine test_infinite_1d

  subroutine test_cuhre(failures)
    integer, intent(inout) :: failures
    type(cubature_result_t) :: result
    real(dp) :: lower2(2), upper2(2), lower3(3), upper3(3)
    real(dp) :: lower4(4), upper4(4), inf

    lower2 = 0.0_dp
    upper2 = 1.0_dp
    call integrate_nd(product_function, lower2, upper2, result, &
                      maxeval=10000, eps_abs=1.0e-10_dp, eps_rel=1.0e-10_dp)
    call check_close('Cuhre 2D product', result%value, 0.25_dp, &
                     1.0e-12_dp, failures)
    call check_true('Cuhre 2D status', result%error_code == 0, failures)

    lower3 = 0.0_dp
    upper3 = 1.0_dp
    call integrate_nd(product_function, lower3, upper3, result, &
                      maxeval=10000, eps_abs=1.0e-10_dp, eps_rel=1.0e-10_dp)
    call check_close('Cuhre 3D product', result%value, 0.125_dp, &
                     1.0e-12_dp, failures)

    lower4 = 0.0_dp
    upper4 = 1.0_dp
    call integrate_nd(product_function, lower4, upper4, result, &
                      maxeval=10000, eps_abs=1.0e-10_dp, eps_rel=1.0e-10_dp)
    call check_close('Cuhre 4D product', result%value, 0.0625_dp, &
                     1.0e-12_dp, failures)

    lower2 = -1.0_dp
    upper2 = 1.0_dp
    call integrate_nd(bivariate_normal, lower2, upper2, result, &
                      maxeval=30000, eps_abs=1.0e-8_dp, eps_rel=1.0e-8_dp)
    call check_close('bivariate normal rectangle', result%value, &
                     0.49797177783920804_dp, 5.0e-8_dp, failures)

    inf = ieee_value(1.0_dp, ieee_positive_inf)
    lower2 = [0.0_dp, -inf]
    upper2 = [inf, inf]
    call integrate_nd(mixed_infinite_density, lower2, upper2, result, &
                      maxeval=50000, eps_abs=1.0e-7_dp, eps_rel=1.0e-7_dp)
    call check_close('mixed infinite Cuhre', result%value, 1.0_dp, &
                     2.0e-6_dp, failures)

    lower2 = [1.0_dp, 0.0_dp]
    upper2 = [0.0_dp, 1.0_dp]
    call integrate_nd(product_function, lower2, upper2, result)
    call check_true('invalid multidimensional bounds', &
                    result%error_code == -2, failures)
  end subroutine test_cuhre

  subroutine test_optimization(failures)
    integer, intent(inout) :: failures
    type(optimization_result_t) :: result
    type(shift_data_t) :: data
    real(dp) :: x(2), lower(2), upper(2), scalar_x(1)

    x = [-1.2_dp, 1.0_dp]
    call optim_lbfgs(rosenbrock, x, result, maxit=500, &
                     eps_f=1.0e-10_dp, eps_g=1.0e-7_dp)
    call check_close('L-BFGS x1', x(1), 1.0_dp, 2.0e-4_dp, failures)
    call check_close('L-BFGS x2', x(2), 1.0_dp, 2.0e-4_dp, failures)
    call check_true('L-BFGS status', result%converged, failures)

    x = [-1.2_dp, 1.0_dp]
    lower = [-2.0_dp, -1.0_dp]
    upper = [2.0_dp, 3.0_dp]
    call optim_lbfgsb(rosenbrock, x, lower, upper, result, maxit=500, &
                      eps_f=1.0e-10_dp, eps_g=1.0e-7_dp)
    call check_close('L-BFGS-B x1', x(1), 1.0_dp, 2.0e-4_dp, failures)
    call check_close('L-BFGS-B x2', x(2), 1.0_dp, 2.0e-4_dp, failures)
    call check_true('L-BFGS-B status', result%converged, failures)

    data%target = 3.25_dp
    scalar_x = 0.0_dp
    call optim_lbfgs(shifted_quadratic, scalar_x, result, &
                     user_data=data, eps_g=1.0e-10_dp)
    call check_close('optimizer user data', scalar_x(1), data%target, &
                     1.0e-7_dp, failures)

    scalar_x = 0.0_dp
    call optim_lbfgsb(bound_quadratic, scalar_x, [-1.0_dp], [2.0_dp], &
                      result, eps_g=1.0e-9_dp)
    call check_close('active upper bound', scalar_x(1), 2.0_dp, &
                     1.0e-7_dp, failures)
  end subroutine test_optimization

  subroutine test_fast_lr(failures)
    integer, intent(inout) :: failures
    type(logistic_fit_t) :: fit
    real(dp) :: x(10,2), y(10)

    x(:,1) = 1.0_dp
    x(:,2) = [-2.0_dp, -1.5_dp, -1.0_dp, -0.5_dp, 0.0_dp, &
              0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp, 2.5_dp]
    y = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
         1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    call fast_lr(x, y, fit, eps_f=1.0e-10_dp, eps_g=1.0e-8_dp)
    call check_close('fastLR intercept', fit%coefficients(1), &
                     -0.27184836_dp, 8.0e-5_dp, failures)
    call check_close('fastLR slope', fit%coefficients(2), &
                     1.08739343_dp, 8.0e-5_dp, failures)
    call check_close('fastLR log likelihood', fit%log_likelihood, &
                     -4.941579983434302_dp, 2.0e-8_dp, failures)
    call check_true('fastLR status', fit%converged, failures)
  end subroutine test_fast_lr

  function power_six(x, user_data) result(value)
    real(dp), intent(in) :: x
    class(*), intent(inout), optional :: user_data
    real(dp) :: value
    value = x**6
  end function power_six

  function beta_pdf(x, user_data) result(value)
    real(dp), intent(in) :: x
    class(*), intent(inout), optional :: user_data
    real(dp) :: value
    value = 660.0_dp*x*x*(1.0_dp - x)**9
  end function beta_pdf

  function gaussian_kernel(x, user_data) result(value)
    real(dp), intent(in) :: x
    class(*), intent(inout), optional :: user_data
    real(dp) :: value
    value = exp(-x*x)
  end function gaussian_kernel

  function exponential_tail(x, user_data) result(value)
    real(dp), intent(in) :: x
    class(*), intent(inout), optional :: user_data
    real(dp) :: value
    value = exp(-x)
  end function exponential_tail

  function product_function(x, user_data) result(value)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: value
    value = product(x)
  end function product_function

  function bivariate_normal(x, user_data) result(value)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: value
    real(dp), parameter :: rho = 0.5_dp
    real(dp) :: quadratic
    quadratic = x(1)**2 - 2.0_dp*rho*x(1)*x(2) + x(2)**2
    value = exp(-quadratic/(2.0_dp*(1.0_dp - rho*rho)))/ &
            (2.0_dp*acos(-1.0_dp)*sqrt(1.0_dp - rho*rho))
  end function bivariate_normal

  function mixed_infinite_density(x, user_data) result(value)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: value
    value = exp(-x(1))*exp(-x(2)*x(2))/sqrt(acos(-1.0_dp))
  end function mixed_infinite_density

  subroutine rosenbrock(x, value, gradient, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value, gradient(:)
    class(*), intent(inout), optional :: user_data
    value = 100.0_dp*(x(2) - x(1)*x(1))**2 + (1.0_dp - x(1))**2
    gradient(1) = -400.0_dp*x(1)*(x(2) - x(1)*x(1)) - &
                  2.0_dp*(1.0_dp - x(1))
    gradient(2) = 200.0_dp*(x(2) - x(1)*x(1))
  end subroutine rosenbrock

  subroutine shifted_quadratic(x, value, gradient, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value, gradient(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: target
    target = 0.0_dp
    if (present(user_data)) then
      select type (user_data)
      type is (shift_data_t)
        target = user_data%target
      end select
    end if
    value = 0.5_dp*(x(1) - target)**2
    gradient(1) = x(1) - target
  end subroutine shifted_quadratic

  subroutine bound_quadratic(x, value, gradient, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value, gradient(:)
    class(*), intent(inout), optional :: user_data
    value = 0.5_dp*(x(1) - 3.0_dp)**2
    gradient(1) = x(1) - 3.0_dp
  end subroutine bound_quadratic

  subroutine check_close(name, actual, expected, tolerance, failures)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: actual, expected, tolerance
    integer, intent(inout) :: failures
    if (abs(actual - expected) > tolerance) then
      write(*,'(a,2(1x,es24.16),a,es12.4)') trim(name), actual, expected, &
                                                  ' tolerance=', tolerance
      failures = failures + 1
    end if
  end subroutine check_close

  subroutine check_true(name, condition, failures)
    character(len=*), intent(in) :: name
    logical, intent(in) :: condition
    integer, intent(inout) :: failures
    if (.not. condition) then
      write(*,'(a)') trim(name)//' failed'
      failures = failures + 1
    end if
  end subroutine check_true

end program test_rcppnumerical
