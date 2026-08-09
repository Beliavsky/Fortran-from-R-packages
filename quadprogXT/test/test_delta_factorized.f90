program test_delta_factorized
  use quadprog_kinds, only: dp
  use quadprog, only: qp_result
  use quadprogxt, only: solve_qp_xt
  use qpxt_test_support, only: assert_true, assert_vector_close
  implicit none

  real(dp) :: dmat(2,2), dvec(2), b0(2)
  real(dp) :: adelta(4,1), bdelta(1)
  type(qp_result) :: a, b

  dmat = 0.0_dp
  dmat(1,1) = 1.0_dp
  dmat(2,2) = 1.0_dp
  dvec = [2.0_dp, -2.0_dp]
  b0 = [0.2_dp, 0.2_dp]
  adelta(:,1) = -1.0_dp
  bdelta(1) = -0.5_dp

  a = solve_qp_xt(dmat, dvec, b0=b0, amat_posneg_delta=adelta, &
    bvec_posneg_delta=bdelta)
  call assert_true(a%succeeded(), 'delta constrained solve')
  call assert_true(sum(abs(a%solution(1:2) - b0)) <= &
    0.5_dp + 2.0e-9_dp, 'delta L1 constraint')

  ! For D=I, R^{-1}=I, so the factorized and ordinary inputs coincide.
  b = solve_qp_xt(dmat, dvec, factorized=.true., b0=b0, &
    amat_posneg_delta=adelta, bvec_posneg_delta=bdelta)
  call assert_true(b%succeeded(), 'factorized delta solve')
  call assert_vector_close(a%solution, b%solution, 2.0e-8_dp, &
    'factorized expanded QP consistency')

  write(*, '(a)') 'PASS test_delta_factorized'
end program test_delta_factorized
