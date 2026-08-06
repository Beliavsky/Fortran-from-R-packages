module test_nonlinear_inequality_problem
  use quadprog_kinds, only: dp
  implicit none
contains
  function objective(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = (1.0_dp - x(1))**2 + 100.0_dp * (x(2) - x(1)**2)**2
  end function objective

  subroutine constraints(x, equality, inequality)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: equality(:), inequality(:)
    allocate(equality(0), inequality(2))
    inequality(1) = x(1)**2 + x(2)**2 - 1.5_dp
    inequality(2) = -x(1) - x(2) - 0.2_dp
  end subroutine constraints
end module test_nonlinear_inequality_problem

program test_nonlinear_inequality
  use quadprog_kinds, only: dp
  use nlcoptim, only: nlc_result, nlc_options, solnl
  use test_nonlinear_inequality_problem, only: objective, constraints
  implicit none
  type(nlc_result) :: fit
  type(nlc_options) :: opt

  opt%central_differences = .true.
  opt%max_iterations = 500
  opt%max_line_search = 45
  call solnl([-1.9_dp, 2.0_dp], objective, fit, confun=constraints, options=opt)
  if (.not. fit%succeeded()) then
    print *, fit%message, fit%x, fit%objective, fit%max_constraint_violation
    error stop 'nonlinear inequality solve failed'
  end if
  if (fit%x(1)**2 + fit%x(2)**2 > 1.50001_dp) error stop 'disk constraint failed'
  if (-fit%x(1) - fit%x(2) > 0.20001_dp) error stop 'halfspace constraint failed'
  if (fit%objective > 0.02_dp) error stop 'poor nonlinear inequality solution'
  print *, 'test_nonlinear_inequality: PASS'
end program test_nonlinear_inequality
