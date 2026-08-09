program test_base_qp
  use quadprog_kinds, only: dp
  use quadprog, only: qp_result, solve_qp
  use quadprogxt, only: solve_qp_xt
  use qpxt_test_support, only: assert_true, assert_vector_close, assert_close
  implicit none

  real(dp) :: dmat(3,3), dvec(3), amat(3,3), bvec(3)
  type(qp_result) :: base, xt

  dmat = 0.0_dp
  dmat(1,1) = 1.0_dp
  dmat(2,2) = 1.0_dp
  dmat(3,3) = 1.0_dp
  dvec = [0.0_dp, 5.0_dp, 0.0_dp]
  amat = reshape([ -4.0_dp, -3.0_dp, 0.0_dp, &
                    2.0_dp,  1.0_dp, 0.0_dp, &
                    0.0_dp, -2.0_dp, 1.0_dp ], [3,3])
  bvec = [-8.0_dp, 2.0_dp, 0.0_dp]

  base = solve_qp(dmat, dvec, amat, bvec)
  xt = solve_qp_xt(dmat, dvec, amat, bvec)
  call assert_true(base%succeeded(), 'base quadprog succeeds')
  call assert_true(xt%succeeded(), 'quadprogXT succeeds')
  call assert_vector_close(xt%solution, base%solution, 2.0e-12_dp, &
    'base QP consistency')
  call assert_close(xt%value, base%value, 2.0e-12_dp, &
    'base QP value consistency')

  write(*, '(a)') 'PASS test_base_qp'
end program test_base_qp
