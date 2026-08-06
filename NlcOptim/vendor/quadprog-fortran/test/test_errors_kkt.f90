! SPDX-License-Identifier: GPL-2.0-or-later
program test_errors_kkt
  use quadprog_kinds, only: dp
  use quadprog, only: qp_result, qp_success, &
    qp_inconsistent_constraints, qp_not_positive_definite, &
    qp_invalid_dimensions, qp_invalid_meq, solve_qp
  use quadprog_test_support, only: check
  implicit none

  real(dp) :: d1(1, 1), dv1(1), a1(1, 2), b1(2)
  real(dp) :: d2(2, 2), dv2(2), a2(2, 4), b2(4), gradient(2)
  real(dp) :: bad_a(3, 1), empty_a(2, 0)
  type(qp_result) :: res

  d1(1, 1) = 1.0_dp
  dv1 = 0.0_dp
  a1 = reshape([1.0_dp, -1.0_dp], [1, 2])
  b1 = [1.0_dp, 0.0_dp]
  res = solve_qp(d1, dv1, a1, b1)
  call check(res%status == qp_inconsistent_constraints, &
    'inconsistent constraints status')

  d2 = reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], [2, 2])
  dv2 = 0.0_dp
  a2 = 0.0_dp
  a2(1, 1) = 1.0_dp
  b2 = -100.0_dp
  res = solve_qp(d2, dv2, a2(:, 1:1), b2(1:1))
  call check(res%status == qp_not_positive_definite, &
    'non-positive-definite status')

  res = solve_qp(d2, dv2, bad_a)
  call check(res%status == qp_invalid_dimensions, 'dimension status')
  res = solve_qp(reshape([1.0_dp], [1, 1]), [0.0_dp], &
    reshape([1.0_dp], [1, 1]), [0.0_dp], meq=2)
  call check(res%status == qp_invalid_meq, 'meq status')

  d2 = reshape([2.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], [2, 2])
  dv2 = [0.2_dp, 0.1_dp]
  a2 = reshape([1.0_dp, 1.0_dp, &
                 1.0_dp, 0.0_dp, &
                 0.0_dp, 1.0_dp, &
                -1.0_dp, 0.0_dp], [2, 4])
  b2 = [1.0_dp, 0.0_dp, 0.0_dp, -0.8_dp]
  res = solve_qp(d2, dv2, a2, b2, meq=1)
  call check(res%status == qp_success, 'KKT problem status')
  call check(abs(sum(res%solution) - 1.0_dp) < 1.0e-11_dp, &
    'KKT equality feasibility')
  call check(minval(matmul(transpose(a2), res%solution) - b2) > &
    -1.0e-11_dp, 'KKT inequality feasibility')
  gradient = matmul(d2, res%solution) - dv2 - &
    matmul(a2, res%lagrangian)
  call check(maxval(abs(gradient)) < 2.0e-10_dp, 'KKT stationarity')

  res = solve_qp(d2, dv2, empty_a)
  call check(res%status == qp_success, 'empty-constraint status')

  write(*, '(a)') 'test_errors_kkt: PASS'
end program test_errors_kkt
