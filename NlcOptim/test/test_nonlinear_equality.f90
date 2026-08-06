module test_nonlinear_equality_problem
  use quadprog_kinds, only: dp
  implicit none
contains
  function objective(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = (x(1) - 1.0_dp)**2 + (x(2) - 2.0_dp)**2
  end function objective

  subroutine constraints(x, equality, inequality)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: equality(:), inequality(:)
    allocate(equality(1), inequality(0))
    equality(1) = x(1)**2 + x(2)**2 - 1.0_dp
  end subroutine constraints

  subroutine jacobian(x, equality_jacobian, inequality_jacobian)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: equality_jacobian(:, :)
    real(dp), allocatable, intent(out) :: inequality_jacobian(:, :)
    allocate(equality_jacobian(1, 2), inequality_jacobian(0, 2))
    equality_jacobian(1, :) = 2.0_dp * x
  end subroutine jacobian
end module test_nonlinear_equality_problem

program test_nonlinear_equality
  use quadprog_kinds, only: dp
  use nlcoptim, only: nlc_result, nlc_options, solnl
  use test_nonlinear_equality_problem, only: objective, constraints, jacobian
  implicit none
  type(nlc_result) :: fit
  type(nlc_options) :: opt
  real(dp) :: expected(2)

  opt%central_differences = .true.
  opt%max_iterations = 300
  call solnl([0.5_dp, 0.5_dp], objective, fit, confun=constraints, options=opt, &
    jacfun=jacobian)
  if (.not. fit%succeeded()) then
    print *, fit%message, fit%max_constraint_violation, fit%kkt_error
    error stop 'nonlinear equality solve failed'
  end if
  expected = [1.0_dp, 2.0_dp] / sqrt(5.0_dp)
  if (maxval(abs(fit%x - expected)) > 5.0e-4_dp) &
    error stop 'wrong nonlinear equality solution'
  if (abs(sum(fit%x**2) - 1.0_dp) > 1.0e-6_dp) &
    error stop 'circle equality not satisfied'
  print *, 'test_nonlinear_equality: PASS'
end program test_nonlinear_equality
