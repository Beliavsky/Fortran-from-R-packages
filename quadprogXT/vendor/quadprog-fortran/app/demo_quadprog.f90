! SPDX-License-Identifier: GPL-2.0-or-later
program demo_quadprog
  use quadprog_kinds, only: dp
  use quadprog, only: qp_result, solve_qp
  implicit none

  real(dp) :: dmat(2, 2), dvec(2), amat(2, 3), bvec(3)
  type(qp_result) :: fit

  dmat = reshape([4.0_dp, -2.0_dp, -2.0_dp, 4.0_dp], [2, 2])
  dvec = [-6.0_dp, 0.0_dp]
  amat = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
    1.0_dp, 1.0_dp], [2, 3])
  bvec = [0.0_dp, 0.0_dp, 2.0_dp]

  fit = solve_qp(dmat, dvec, amat, bvec)
  write(*, '(a,i0,2a)') 'status: ', fit%status, ' (', fit%message // ')'
  write(*, '(a,*(f10.5,1x))') 'solution: ', fit%solution
  write(*, '(a,*(f10.5,1x))') 'multipliers: ', fit%lagrangian
  write(*, '(a,f10.5)') 'objective: ', fit%value
end program demo_quadprog
