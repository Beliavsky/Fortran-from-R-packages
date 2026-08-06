module test_unconstrained_problem
  use quadprog_kinds, only: dp
  implicit none
contains
  function objective(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = (x(1) - 2.0_dp)**2 + 3.0_dp * (x(2) + 1.0_dp)**2
  end function objective

  subroutine gradient(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    g = [2.0_dp * (x(1) - 2.0_dp), 6.0_dp * (x(2) + 1.0_dp)]
  end subroutine gradient
end module test_unconstrained_problem

program test_unconstrained
  use quadprog_kinds, only: dp
  use nlcoptim, only: nlc_result, solnl
  use test_unconstrained_problem, only: objective, gradient
  implicit none
  type(nlc_result) :: fit

  call solnl([8.0_dp, -7.0_dp], objective, fit, gradfun=gradient)
  if (.not. fit%succeeded()) error stop 'unconstrained solve failed'
  if (maxval(abs(fit%x - [2.0_dp, -1.0_dp])) > 1.0e-5_dp) &
    error stop 'wrong unconstrained solution'
  if (fit%objective > 1.0e-10_dp) error stop 'wrong unconstrained objective'
  print *, 'test_unconstrained: PASS'
end program test_unconstrained
