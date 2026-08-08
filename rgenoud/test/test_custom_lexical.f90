program test_custom_lexical
  use rgenoud, only : dp, genoud_options, genoud_result, genoud_optimize_lexical
  implicit none
  type(genoud_options) :: opt
  type(genoud_result) :: res
  real(dp) :: lower(1), upper(1)

  lower = -3.0_dp
  upper = 3.0_dp
  opt%pop_size = 60
  opt%max_generations = 60
  opt%wait_generations = 8
  opt%use_bfgs = .false.
  opt%gradient_check = .false.
  opt%seed = 822
  call genoud_optimize_lexical(obj, 2, lower, upper, opt, res, comparator=compare_second_first)
  if (res%fit(2) > 1.0e-3_dp) error stop "custom lexical comparator failed"
contains
  subroutine obj(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)
    f(1) = (x(1) - 2.0_dp)**2
    f(2) = (x(1) + 1.0_dp)**2
  end subroutine obj
  logical function compare_second_first(a, b) result(better)
    real(dp), intent(in) :: a(:), b(:)
    if (abs(a(2) - b(2)) > 1.0e-14_dp) then
      better = a(2) < b(2)
    else
      better = a(1) < b(1)
    end if
  end function compare_second_first
end program test_custom_lexical
