! SPDX-License-Identifier: BSD-2-Clause
program example_linkage_methods
  use fastcluster
  implicit none

  real(dp), allocatable :: distances(:, :)
  real(dp) :: x(5, 2)
  type(hclust_result) :: result
  character(len=10) :: methods(8)
  integer :: i

  x = reshape([0.0_dp, 1.0_dp, 0.0_dp, 5.0_dp, 6.0_dp, &
               0.0_dp, 0.0_dp, 2.0_dp, 5.0_dp, 4.0_dp], shape(x))
  call pairwise_distances(x, 'euclidean', distances)
  methods = [character(len=10) :: 'single', 'complete', 'average', 'mcquitty', &
    'ward.D', 'centroid', 'median', 'ward.D2']

  do i = 1, size(methods)
    call hclust(distances, trim(methods(i)), result)
    if (.not. result%ok()) error stop result%message
    write (*, '(a10,2x,*(f10.5,1x))') trim(methods(i)), result%height
  end do
end program example_linkage_methods
