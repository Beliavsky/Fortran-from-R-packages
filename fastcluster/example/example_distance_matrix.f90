! SPDX-License-Identifier: BSD-2-Clause
program example_distance_matrix
  use fastcluster
  implicit none

  real(dp) :: d(6)
  type(hclust_result) :: result
  integer :: i

  ! R-compatible condensed order for four observations:
  ! d(2,1), d(3,1), d(4,1), d(3,2), d(4,2), d(4,3).
  d = [1.0_dp, 2.0_dp, sqrt(50.0_dp), sqrt(5.0_dp), sqrt(41.0_dp), sqrt(34.0_dp)]
  call hclust(d, 4, 'complete', result)
  if (.not. result%ok()) error stop result%message

  write (*, '(a)') 'Complete linkage from condensed distances'
  do i = 1, result%n - 1
    write (*, '(2i6,f13.6)') result%merge(i, :), result%height(i)
  end do
end program example_distance_matrix
