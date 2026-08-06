module basic_nonlinear_problem
  use quadprog_kinds, only: dp
  implicit none
contains
  function objective(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = exp(product(x))
  end function objective

  subroutine constraints(x, equality, inequality)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: equality(:), inequality(:)
    allocate(equality(3), inequality(0))
    equality(1) = sum(x**2) - 10.0_dp
    equality(2) = x(2) * x(3) - 5.0_dp * x(4) * x(5)
    equality(3) = x(1)**3 + x(2)**3 + 1.0_dp
  end subroutine constraints
end module basic_nonlinear_problem

program basic_nonlinear
  use quadprog_kinds, only: dp
  use nlcoptim, only: nlc_result, nlc_options, solnl
  use basic_nonlinear_problem, only: objective, constraints
  implicit none
  type(nlc_result) :: fit
  type(nlc_options) :: opt

  opt%central_differences = .true.
  opt%max_iterations = 600
  call solnl([-2.0_dp, 2.0_dp, 2.0_dp, -1.0_dp, -1.0_dp], objective, fit, &
    confun=constraints, options=opt)
  print '(a,l1)', 'success: ', fit%succeeded()
  print '(a,es16.8)', 'objective: ', fit%objective
  print '(a,*(f12.6,1x))', 'x: ', fit%x
  print '(a,es12.4)', 'max constraint violation: ', fit%max_constraint_violation
  print '(a,i0)', 'iterations: ', fit%iterations
end program basic_nonlinear
