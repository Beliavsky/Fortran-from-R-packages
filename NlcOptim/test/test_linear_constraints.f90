module test_linear_problem
  use quadprog_kinds, only: dp
  implicit none
contains
  function objective(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = (x(1) - 0.8_dp)**2 + (x(2) - 0.2_dp)**2
  end function objective
end module test_linear_problem

program test_linear_constraints
  use quadprog_kinds, only: dp
  use nlcoptim, only: nlc_result, solnl
  use test_linear_problem, only: objective
  implicit none
  type(nlc_result) :: fit
  real(dp) :: aeq(1, 2), beq(1), lb(2), ub(2)

  aeq = reshape([1.0_dp, 1.0_dp], [1, 2])
  beq = [1.0_dp]
  lb = [0.0_dp, 0.0_dp]
  ub = [0.6_dp, 1.0_dp]
  call solnl([0.1_dp, 0.9_dp], objective, fit, aeq=aeq, beq=beq, lb=lb, ub=ub)
  if (.not. fit%succeeded()) then
    print *, fit%message
    error stop 'linear constrained solve failed'
  end if
  if (maxval(abs(fit%x - [0.6_dp, 0.4_dp])) > 2.0e-4_dp) &
    error stop 'wrong linearly constrained solution'
  if (abs(sum(fit%x) - 1.0_dp) > 1.0e-7_dp) error stop 'equality not satisfied'
  if (fit%x(1) > 0.600001_dp) error stop 'upper bound not satisfied'
  print *, 'test_linear_constraints: PASS'
end program test_linear_constraints
