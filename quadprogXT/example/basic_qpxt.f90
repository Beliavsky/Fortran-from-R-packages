program basic_qpxt
  use quadprog_kinds, only: dp
  use quadprog, only: qp_result
  use quadprogxt, only: solve_qp_xt
  implicit none

  real(dp) :: dmat(3,3), dvec(3), amat(3,3), bvec(3)
  type(qp_result) :: fit

  dmat = 0.0_dp
  dmat(1,1) = 1.0_dp
  dmat(2,2) = 1.0_dp
  dmat(3,3) = 1.0_dp
  dvec = [0.0_dp, 5.0_dp, 0.0_dp]
  amat = reshape([ -4.0_dp, -3.0_dp, 0.0_dp, &
                    2.0_dp,  1.0_dp, 0.0_dp, &
                    0.0_dp, -2.0_dp, 1.0_dp ], [3,3])
  bvec = [-8.0_dp, 2.0_dp, 0.0_dp]

  fit = solve_qp_xt(dmat, dvec, amat, bvec)
  if (.not. fit%succeeded()) error stop fit%message
  write(*, '(a,3f12.6)') 'solution:', fit%solution(1:3)
  write(*, '(a,es14.6)') 'objective:', fit%value
end program basic_qpxt
