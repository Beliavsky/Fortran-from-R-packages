! SPDX-License-Identifier: GPL-2.0-or-later
program compact_qp
  use quadprog_kinds, only: dp
  use quadprog, only: qp_result, solve_qp_compact
  implicit none

  real(dp) :: dmat(3, 3), dvec(3), values(2, 3), bvec(3)
  integer :: indices(3, 3)
  type(qp_result) :: fit

  dmat = 0.0_dp
  dmat(1, 1) = 1.0_dp
  dmat(2, 2) = 1.0_dp
  dmat(3, 3) = 1.0_dp
  dvec = [0.0_dp, 5.0_dp, 0.0_dp]

  ! Column j stores its nonzero values in values(:,j).  indices(1,j)
  ! is the count, followed by the corresponding variable indices.
  values = reshape([-4.0_dp, -3.0_dp, &
                     2.0_dp,  1.0_dp, &
                    -2.0_dp,  1.0_dp], [2, 3])
  indices = reshape([2, 1, 2, &
                     2, 1, 2, &
                     2, 2, 3], [3, 3])
  bvec = [-8.0_dp, 2.0_dp, 0.0_dp]

  fit = solve_qp_compact(dmat, dvec, values, indices, bvec)
  if (.not. fit%succeeded()) error stop fit%message

  write(*, '(a,*(f10.6,1x))') 'solution: ', fit%solution
  write(*, '(a,f10.6)') 'objective: ', fit%value
end program compact_qp
