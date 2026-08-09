program test_full_problem
  use quadprog_kinds, only: dp
  use quadprog, only: qp_result
  use quadprogxt, only: qpxt_problem, build_qp_xt, solve_qp_xt
  use qpxt_test_support, only: assert_true, assert_vector_close
  implicit none

  integer, parameter :: n = 3
  real(dp) :: dmat(n,n), dvec(n), amat(n,2*n), bvec(2*n)
  real(dp) :: apos(2*n,1), adelta(2*n,1), b0(n)
  real(dp) :: bp(1), bd(1), dpn(2*n), dd(2*n)
  type(qp_result) :: dense, compact
  type(qpxt_problem) :: p
  integer :: i

  dmat = 0.0_dp
  do i = 1, n
    dmat(i,i) = 1.0_dp + 0.2_dp * real(i - 1, dp)
  end do
  dvec = [0.7_dp, -0.3_dp, 0.5_dp]
  amat = 0.0_dp
  do i = 1, n
    amat(i,i) = 1.0_dp
    amat(i,n+i) = -1.0_dp
  end do
  bvec = -2.0_dp
  apos(:,1) = -1.0_dp
  adelta(:,1) = -1.0_dp
  bp = -1.4_dp
  bd = -0.7_dp
  b0 = [0.1_dp, -0.1_dp, 0.2_dp]
  dpn = -0.03_dp
  dd = -0.02_dp

  p = build_qp_xt(dmat, dvec, amat, bvec, amat_posneg=apos, &
    bvec_posneg=bp, dvec_posneg=dpn, b0=b0, &
    amat_posneg_delta=adelta, bvec_posneg_delta=bd, &
    dvec_posneg_delta=dd)
  call assert_true(p%succeeded(), 'full build succeeds')
  call assert_true(size(p%dmat,1) == 3*n, 'full problem has 3n variables')
  call assert_true(p%compact, 'compact is default')

  compact = solve_qp_xt(dmat, dvec, amat, bvec, &
    amat_posneg=apos, bvec_posneg=bp, dvec_posneg=dpn, b0=b0, &
    amat_posneg_delta=adelta, bvec_posneg_delta=bd, &
    dvec_posneg_delta=dd, compact=.true.)
  dense = solve_qp_xt(dmat, dvec, amat, bvec, &
    amat_posneg=apos, bvec_posneg=bp, dvec_posneg=dpn, b0=b0, &
    amat_posneg_delta=adelta, bvec_posneg_delta=bd, &
    dvec_posneg_delta=dd, compact=.false.)
  call assert_true(compact%succeeded(), 'full compact solve')
  call assert_true(dense%succeeded(), 'full dense solve')
  call assert_vector_close(compact%solution, dense%solution, 3.0e-8_dp, &
    'dense and compact full problem agree')

  write(*, '(a)') 'PASS test_full_problem'
end program test_full_problem
