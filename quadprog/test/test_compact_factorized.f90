! SPDX-License-Identifier: GPL-2.0-or-later
program test_compact_factorized
  use quadprog_kinds, only: dp
  use quadprog, only: qp_result, qp_success, qp_invalid_compact_index, &
    solve_qp, solve_qp_compact
  use quadprog_test_support, only: check, check_close, check_vector_close
  implicit none

  real(dp) :: d(2, 2), rinv(2, 2), dvec(2), ac(2, 3), b(3)
  real(dp) :: afull(2, 3), empty_a(2, 0)
  integer :: aind(3, 3), bad_aind(3, 3)
  type(qp_result) :: compact_res, factor_res, unconstrained_res

  d = reshape([4.0_dp, -2.0_dp, -2.0_dp, 4.0_dp], [2, 2])
  rinv = 0.0_dp
  rinv(1, 1) = 0.5_dp
  rinv(1, 2) = 0.5_dp / sqrt(3.0_dp)
  rinv(2, 2) = 1.0_dp / sqrt(3.0_dp)
  dvec = [-6.0_dp, 0.0_dp]
  ac = reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
    [2, 3])
  aind = reshape([1, 1, 0, 1, 2, 0, 2, 1, 2], [3, 3])
  b = [0.0_dp, 0.0_dp, 2.0_dp]

  compact_res = solve_qp_compact(d, dvec, ac, aind, b)
  call check(compact_res%status == qp_success, 'compact status')
  call check_vector_close(compact_res%solution, [0.5_dp, 1.5_dp], &
    1.0e-12_dp, 'compact solution')
  call check_close(compact_res%value, 6.5_dp, 1.0e-12_dp, &
    'compact value')

  factor_res = solve_qp_compact(rinv, dvec, ac, aind, b, &
    factorized=.true.)
  call check(factor_res%status == qp_success, 'compact factorized status')
  call check_vector_close(factor_res%solution, compact_res%solution, &
    2.0e-12_dp, 'compact factorized solution')

  afull = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
    1.0_dp, 1.0_dp], [2, 3])
  factor_res = solve_qp(rinv, dvec, afull, b, factorized=.true.)
  call check(factor_res%status == qp_success, 'factorized status')
  call check_vector_close(factor_res%solution, compact_res%solution, &
    2.0e-12_dp, 'factorized solution')
  call check_vector_close(factor_res%unconstrained_solution, &
    [-2.0_dp, -1.0_dp], 2.0e-12_dp, 'factorized unconstrained')

  unconstrained_res = solve_qp(d, dvec, empty_a)
  call check(unconstrained_res%succeeded(), 'unconstrained status')
  call check_vector_close(unconstrained_res%solution, &
    [-2.0_dp, -1.0_dp], 2.0e-12_dp, 'unconstrained result')
  call check_close(unconstrained_res%value, -6.0_dp, 2.0e-12_dp, &
    'unconstrained value')
  unconstrained_res = solve_qp(rinv, dvec, empty_a, factorized=.true.)
  call check_vector_close(unconstrained_res%solution, &
    [-2.0_dp, -1.0_dp], 2.0e-12_dp, 'factorized unconstrained result')

  bad_aind = aind
  bad_aind(2, 1) = 3
  compact_res = solve_qp_compact(d, dvec, ac, bad_aind, b)
  call check(compact_res%status == qp_invalid_compact_index, &
    'invalid compact index status')

  write(*, '(a)') 'test_compact_factorized: PASS'
end program test_compact_factorized
