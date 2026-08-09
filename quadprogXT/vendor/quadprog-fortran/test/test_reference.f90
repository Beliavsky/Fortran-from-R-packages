! SPDX-License-Identifier: GPL-2.0-or-later
program test_reference
  use quadprog_kinds, only: dp
  use quadprog, only: qp_result, qp_success, solve_qp
  use quadprog_test_support, only: check, check_close, check_vector_close
  implicit none

  real(dp) :: d2(2, 2), dvec2(2), a2(2, 3), b2(3)
  real(dp) :: d3(3, 3), dvec3(3), a3(3, 3), b3(3)
  type(qp_result) :: res

  d2 = reshape([4.0_dp, -2.0_dp, -2.0_dp, 4.0_dp], [2, 2])
  dvec2 = [-6.0_dp, 0.0_dp]
  a2 = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
    1.0_dp, 1.0_dp], [2, 3])
  b2 = [0.0_dp, 0.0_dp, 2.0_dp]

  res = solve_qp(d2, dvec2, a2, b2)
  call check(res%status == qp_success, 'two-variable reference status')
  call check_vector_close(res%solution, [0.5_dp, 1.5_dp], 1.0e-12_dp, &
    'two-variable reference solution')
  call check_vector_close(res%unconstrained_solution, &
    [-2.0_dp, -1.0_dp], 1.0e-12_dp, 'unconstrained solution')
  call check_close(res%value, 6.5_dp, 1.0e-12_dp, 'reference value')
  call check_vector_close(res%lagrangian, [0.0_dp, 0.0_dp, 5.0_dp], &
    1.0e-12_dp, 'reference multipliers')
  call check(all(res%active_set == [3]), 'reference active set')
  call check(all(res%iterations == [2, 0]), 'reference iterations')

  d3 = 0.0_dp
  d3(1, 1) = 1.0_dp
  d3(2, 2) = 1.0_dp
  d3(3, 3) = 1.0_dp
  dvec3 = [0.0_dp, 5.0_dp, 0.0_dp]
  a3 = reshape([-4.0_dp, -3.0_dp, 0.0_dp, &
                 2.0_dp,  1.0_dp, 0.0_dp, &
                 0.0_dp, -2.0_dp, 1.0_dp], [3, 3])
  b3 = [-8.0_dp, 2.0_dp, 0.0_dp]

  res = solve_qp(d3, dvec3, a3, b3)
  call check(res%succeeded(), 'three-variable reference status')
  call check_vector_close(res%solution, &
    [10.0_dp / 21.0_dp, 22.0_dp / 21.0_dp, 44.0_dp / 21.0_dp], &
    2.0e-12_dp, 'three-variable reference solution')
  call check_close(res%value, -50.0_dp / 21.0_dp, 2.0e-12_dp, &
    'three-variable reference value')
  call check_vector_close(res%lagrangian, &
    [0.0_dp, 5.0_dp / 21.0_dp, 44.0_dp / 21.0_dp], &
    2.0e-12_dp, 'three-variable multipliers')
  call check(all(res%active_set == [3, 2]), 'three-variable active set')

  write(*, '(a)') 'test_reference: PASS'
end program test_reference
