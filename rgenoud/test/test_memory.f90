program test_memory
  use rgenoud, only : dp, genoud_options, genoud_result, genoud_optimize
  implicit none
  type(genoud_options) :: opt
  type(genoud_result) :: res
  real(dp) :: lower(1), upper(1), starts(1, 6)

  lower = -1.0_dp
  upper = 1.0_dp
  starts = 0.25_dp
  opt%pop_size = 12
  opt%max_generations = 5
  opt%wait_generations = 2
  opt%memory_matrix = .true.
  opt%gradient_check = .false.
  opt%use_bfgs = .false.
  opt%seed = 42
  call genoud_optimize(obj, lower, upper, opt, res, starting_values=starts)
  if (res%unique_evaluations >= res%evaluations) then
    error stop "memory cache did not eliminate repeated evaluation"
  end if
contains
  real(dp) function obj(x) result(f)
    real(dp), intent(in) :: x(:)
    f = x(1)**2
  end function obj
end program test_memory
