module test_error_problem
  use quadprog_kinds, only: dp
  implicit none
contains
  function objective(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = sum(x**2)
  end function objective
end module test_error_problem

program test_qp_and_errors
  use quadprog_kinds, only: dp
  use quadprog, only: qp_result, solve_qp
  use nlcoptim, only: nlc_result, nlc_invalid_input, solnl
  use test_error_problem, only: objective
  implicit none
  type(qp_result) :: qp
  type(nlc_result) :: fit
  real(dp) :: d(2, 2), dv(2), am(2, 1), bv(1)

  d = reshape([2.0_dp, 0.0_dp, 0.0_dp, 2.0_dp], [2, 2])
  dv = [2.0_dp, 4.0_dp]
  am(:, 1) = [1.0_dp, 1.0_dp]
  bv = [3.0_dp]
  qp = solve_qp(d, dv, am, bv, meq=1)
  if (.not. qp%succeeded()) error stop 'vendored quadprog failed'
  if (maxval(abs(qp%solution - [1.0_dp, 2.0_dp])) > 1.0e-10_dp) &
    error stop 'wrong QP solution'

  call solnl([0.0_dp, 0.0_dp], objective, fit, lb=[1.0_dp], ub=[2.0_dp, 2.0_dp])
  if (fit%status /= nlc_invalid_input) error stop 'invalid dimensions not rejected'
  print *, 'test_qp_and_errors: PASS'
end program test_qp_and_errors
