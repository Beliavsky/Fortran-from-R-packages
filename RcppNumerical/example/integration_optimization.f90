program integration_optimization
  use rcppnumerical
  implicit none

  type(integration_result_t) :: one_d
  type(cubature_result_t) :: multi_d
  type(optimization_result_t) :: opt
  real(dp) :: lower(2), upper(2), x(2)

  call integrate_1d(beta_pdf, 0.3_dp, 0.8_dp, one_d)
  write(*,'(a,f14.10,a,es12.4)') 'Beta probability = ', one_d%value, &
                                 ', estimated error = ', one_d%error_estimate

  lower = -1.0_dp
  upper = 1.0_dp
  call integrate_nd(bivariate_normal, lower, upper, multi_d, maxeval=20000)
  write(*,'(a,f14.10,a,es12.4)') 'Bivariate normal probability = ', &
    multi_d%value, ', estimated error = ', multi_d%error_estimate

  x = [-1.2_dp, 1.0_dp]
  call optim_lbfgs(rosenbrock, x, opt)
  write(*,'(a,2f14.8,a,es12.4)') 'Rosenbrock optimum = ', x, &
                                 ', objective = ', opt%value

contains

  function beta_pdf(x, user_data) result(value)
    real(dp), intent(in) :: x
    class(*), intent(inout), optional :: user_data
    real(dp) :: value
    value = 660.0_dp*x*x*(1.0_dp - x)**9
  end function beta_pdf

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

  subroutine rosenbrock(x, value, gradient, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value, gradient(:)
    class(*), intent(inout), optional :: user_data
    value = 100.0_dp*(x(2) - x(1)*x(1))**2 + (1.0_dp - x(1))**2
    gradient(1) = -400.0_dp*x(1)*(x(2) - x(1)*x(1)) - &
                  2.0_dp*(1.0_dp - x(1))
    gradient(2) = 200.0_dp*(x(2) - x(1)*x(1))
  end subroutine rosenbrock

end program integration_optimization
