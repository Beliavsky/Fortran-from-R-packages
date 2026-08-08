program test_integer
  use rgenoud, only : dp, genoud_options, genoud_result, genoud_optimize
  implicit none
  type(genoud_options) :: opt
  type(genoud_result) :: res
  real(dp) :: lower(2), upper(2)

  lower = [-10.0_dp, -10.0_dp]
  upper = [10.0_dp, 10.0_dp]
  opt%integer_parameters = .true.
  opt%pop_size = 80
  opt%max_generations = 80
  opt%wait_generations = 10
  opt%seed = 287
  call genoud_optimize(integer_objective, lower, upper, opt, res)
  if (res%fit(1) > 0.0_dp) error stop "integer optimization test failed"
  if (any(abs(res%par - [3.0_dp, -2.0_dp]) > 0.0_dp)) then
    error stop "integer parameters test failed"
  end if
contains
  real(dp) function integer_objective(x) result(f)
    real(dp), intent(in) :: x(:)
    f = (x(1) - 3.0_dp)**2 + (x(2) + 2.0_dp)**2
  end function integer_objective
end program test_integer
