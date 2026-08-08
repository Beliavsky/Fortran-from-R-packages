program test_maximize
  use rgenoud, only : dp, genoud_options, genoud_result, genoud_optimize
  implicit none
  type(genoud_options) :: opt
  type(genoud_result) :: res
  real(dp) :: lower(1), upper(1)

  lower = -4.0_dp
  upper = 4.0_dp
  opt%maximize = .true.
  opt%pop_size = 40
  opt%max_generations = 60
  opt%wait_generations = 8
  opt%solution_tolerance = 1.0e-7_dp
  opt%seed = 91
  call genoud_optimize(sine_objective, lower, upper, opt, res)
  if (res%fit(1) < 0.999999_dp) error stop "maximization test failed"
contains
  real(dp) function sine_objective(x) result(f)
    real(dp), intent(in) :: x(:)
    f = sin(x(1))
  end function sine_objective
end program test_maximize
