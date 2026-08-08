program rosenbrock_example
  use rgenoud, only : dp, genoud_options, genoud_result, genoud_optimize
  implicit none
  type(genoud_options) :: opt
  type(genoud_result) :: res
  real(dp) :: lower(2), upper(2)
  lower = -5.0_dp
  upper = 5.0_dp
  opt%pop_size = 80
  opt%max_generations = 80
  opt%wait_generations = 10
  opt%boundary_enforcement = 2
  opt%seed = 1234
  call genoud_optimize(rosen, lower, upper, opt, res)
  write(*, '(a,es15.7)') 'fitness: ', res%fit(1)
  write(*, '(a,*(1x,es15.7))') 'parameters:', res%par
contains
  real(dp) function rosen(x) result(f)
    real(dp), intent(in) :: x(:)
    f = 100.0_dp * (x(2) - x(1)**2)**2 + (1.0_dp - x(1))**2
  end function rosen
end program rosenbrock_example
