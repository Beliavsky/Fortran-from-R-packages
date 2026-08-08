program test_lexical
  use rgenoud, only : dp, genoud_options, genoud_result, genoud_optimize_lexical
  implicit none
  type(genoud_options) :: opt
  type(genoud_result) :: res
  real(dp) :: lower(2), upper(2)

  lower = -4.0_dp
  upper = 4.0_dp
  opt%pop_size = 100
  opt%max_generations = 100
  opt%wait_generations = 15
  opt%gradient_check = .false.
  opt%use_bfgs = .false.
  opt%seed = 531
  call genoud_optimize_lexical(lexobj, 2, lower, upper, opt, res)
  if (res%fit(1) > 5.0e-3_dp) error stop "lexical primary criterion failed"
  if (abs(res%par(1) - 1.0_dp) > 0.08_dp) error stop "lexical x1 failed"
contains
  subroutine lexobj(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(1) = (x(1) - 1.0_dp)**2
    f(2) = (x(2) + 2.0_dp)**2
  end subroutine lexobj
end program test_lexical
