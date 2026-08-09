program test_absolute_value
  use quadprog_kinds, only: dp
  use quadprog, only: qp_result
  use quadprogxt, only: solve_qp_xt
  use qpxt_test_support, only: assert_true, assert_close
  implicit none

  real(dp) :: dmat(2,2), dvec(2)
  real(dp) :: abs_a(4,1), babs(1), l1(4)
  type(qp_result) :: uncon, constrained, penalized

  dmat = 0.0_dp
  dmat(1,1) = 1.0_dp
  dmat(2,2) = 1.0_dp
  dvec = [2.0_dp, 1.0_dp]

  uncon = solve_qp_xt(dmat, dvec, compact=.false.)
  call assert_true(uncon%succeeded(), 'unconstrained solve')
  call assert_close(uncon%solution(1), 2.0_dp, 1.0e-12_dp, 'uncon x1')
  call assert_close(uncon%solution(2), 1.0_dp, 1.0e-12_dp, 'uncon x2')

  abs_a(:,1) = -1.0_dp
  babs(1) = -1.0_dp
  constrained = solve_qp_xt(dmat, dvec, amat_posneg=abs_a, &
    bvec_posneg=babs)
  call assert_true(constrained%succeeded(), 'absolute constraint solve')
  call assert_true(sum(abs(constrained%solution(1:2))) <= &
    1.0_dp + 2.0e-9_dp, 'L1 ball is enforced')
  call assert_close(constrained%solution(1), 1.0_dp, 2.0e-7_dp, &
    'L1 constrained x1')
  call assert_close(constrained%solution(2), 0.0_dp, 2.0e-7_dp, &
    'L1 constrained x2')

  l1 = -0.5_dp
  penalized = solve_qp_xt(dmat, dvec, dvec_posneg=l1)
  call assert_true(penalized%succeeded(), 'L1 penalty solve')
  call assert_true(sum(abs(penalized%solution(1:2))) < &
    sum(abs(uncon%solution(1:2))), 'L1 penalty shrinks solution')
  call assert_close(penalized%solution(1), 1.5_dp, 2.0e-7_dp, &
    'L1 penalized x1')
  call assert_close(penalized%solution(2), 0.5_dp, 2.0e-7_dp, &
    'L1 penalized x2')

  write(*, '(a)') 'PASS test_absolute_value'
end program test_absolute_value
