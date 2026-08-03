! SPDX-License-Identifier: GPL-2.0-or-later
program example_hierarchy
  use cluster, only: dp, hierarchy_result, agnes, diana
  implicit none

  real(dp) :: x(6,2)
  type(hierarchy_result) :: result
  integer :: i

  x = reshape([0.0_dp,0.2_dp,-0.1_dp,4.8_dp,5.0_dp,5.2_dp, &
               0.0_dp,-0.1_dp,0.2_dp,5.1_dp,4.9_dp,5.0_dp], [6,2])

  call agnes(x, result, method='average')
  write(*,'(a,f8.4)') 'Agglomerative coefficient: ', result%coefficient
  do i = 1, size(result%height)
    write(*,'(i3,2i6,f12.5)') i, result%merge(i,:), result%height(i)
  end do

  call diana(x, result)
  write(*,'(a,f8.4)') 'Divisive coefficient: ', result%coefficient
  write(*,'(a,*(i0,1x))') 'DIANA order: ', result%order
end program example_hierarchy
