program test_continuous
  use rgenoud, only : dp, genoud_options, genoud_result, genoud_optimize
  implicit none
  type(genoud_options) :: opt
  type(genoud_result) :: res
  real(dp) :: lower(2), upper(2)

  lower = -5.0_dp
  upper = 5.0_dp
  opt%pop_size = 60
  opt%max_generations = 80
  opt%wait_generations = 12
  opt%solution_tolerance = 1.0e-6_dp
  opt%bfgs_gtol = 1.0e-9_dp
  opt%seed = 12345
  opt%boundary_enforcement = 2
  call genoud_optimize(rosenbrock, lower, upper, opt, res)
  if (res%fit(1) > 1.0e-8_dp) error stop "continuous Rosenbrock test failed"
  if (maxval(abs(res%par - [1.0_dp, 1.0_dp])) > 1.0e-3_dp) then
    error stop "continuous parameters test failed"
  end if
contains
  real(dp) function rosenbrock(x) result(f)
    real(dp), intent(in) :: x(:)
    f = 100.0_dp * (x(2) - x(1)**2)**2 + (1.0_dp - x(1))**2
  end function rosenbrock
end program test_continuous
