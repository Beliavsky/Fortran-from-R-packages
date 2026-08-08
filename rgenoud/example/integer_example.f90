program integer_example
  use rgenoud, only : dp, genoud_options, genoud_result, genoud_optimize
  implicit none
  type(genoud_options) :: opt
  type(genoud_result) :: res
  real(dp) :: lower(3), upper(3)
  lower = -10.0_dp
  upper = 10.0_dp
  opt%integer_parameters = .true.
  opt%pop_size = 100
  opt%seed = 777
  call genoud_optimize(target, lower, upper, opt, res)
  write(*, '(a,es15.7)') 'fitness: ', res%fit(1)
  write(*, '(a,*(1x,f6.0))') 'integer parameters:', res%par
contains
  real(dp) function target(x) result(f)
    real(dp), intent(in) :: x(:)
    f = (x(1) - 2.0_dp)**2 + (x(2) + 4.0_dp)**2 + (x(3) - 7.0_dp)**2
  end function target
end program integer_example
