! SPDX-License-Identifier: GPL-2.0-or-later
program portfolio_qp
  use quadprog_kinds, only: dp
  use quadprog, only: qp_result, solve_qp
  implicit none

  integer, parameter :: n = 4
  real(dp) :: covariance(n, n), expected_return(n)
  real(dp) :: constraints(n, n + 1), bounds(n + 1)
  type(qp_result) :: fit
  integer :: i

  covariance = reshape([ &
    0.040_dp, 0.012_dp, 0.008_dp, 0.006_dp, &
    0.012_dp, 0.090_dp, 0.015_dp, 0.010_dp, &
    0.008_dp, 0.015_dp, 0.160_dp, 0.020_dp, &
    0.006_dp, 0.010_dp, 0.020_dp, 0.062_dp], [n, n])
  expected_return = [0.050_dp, 0.075_dp, 0.110_dp, 0.065_dp]

  constraints = 0.0_dp
  constraints(:, 1) = 1.0_dp
  do i = 1, n
    constraints(i, i + 1) = 1.0_dp
  end do
  bounds = 0.0_dp
  bounds(1) = 1.0_dp

  ! Minimize 0.5*w'Cov*w - expected_return'w, subject to sum(w)=1, w>=0.
  fit = solve_qp(covariance, expected_return, constraints, bounds, meq=1)
  if (.not. fit%succeeded()) error stop fit%message

  write(*, '(a,*(f9.5,1x))') 'weights: ', fit%solution
  write(*, '(a,f10.6)') 'portfolio return: ', &
    dot_product(expected_return, fit%solution)
  write(*, '(a,f10.6)') 'portfolio variance: ', &
    dot_product(fit%solution, matmul(covariance, fit%solution))
end program portfolio_qp
