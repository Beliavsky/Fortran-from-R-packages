program l1_qpxt
  use quadprog_kinds, only: dp
  use quadprog, only: qp_result
  use quadprogxt, only: solve_qp_xt
  implicit none

  real(dp) :: dmat(2,2), dvec(2), aabs(4,1), babs(1), penalty(4)
  type(qp_result) :: constrained, penalized

  dmat = 0.0_dp
  dmat(1,1) = 1.0_dp
  dmat(2,2) = 1.0_dp
  dvec = [2.0_dp, 1.0_dp]
  aabs(:,1) = -1.0_dp
  babs = -1.0_dp
  penalty = -0.5_dp

  constrained = solve_qp_xt(dmat, dvec, amat_posneg=aabs, &
    bvec_posneg=babs)
  penalized = solve_qp_xt(dmat, dvec, dvec_posneg=penalty)

  if (.not. constrained%succeeded()) error stop constrained%message
  if (.not. penalized%succeeded()) error stop penalized%message
  write(*, '(a,2f12.6)') 'L1 constrained b:', constrained%solution(1:2)
  write(*, '(a,f12.6)') 'sum(abs(b)):', sum(abs(constrained%solution(1:2)))
  write(*, '(a,2f12.6)') 'L1 penalized b:', penalized%solution(1:2)
end program l1_qpxt
