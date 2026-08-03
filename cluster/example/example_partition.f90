! SPDX-License-Identifier: GPL-2.0-or-later
program example_partition
  use cluster, only: dp, partition_result, pam, clara, fanny
  implicit none

  real(dp) :: x(8,2)
  type(partition_result) :: result
  integer :: i

  x = reshape([0.0_dp,0.2_dp,-0.1_dp,0.1_dp,5.0_dp,5.2_dp,4.9_dp,5.1_dp, &
               0.0_dp,-0.1_dp,0.2_dp,0.1_dp,5.0_dp,4.8_dp,5.2_dp,5.1_dp], [8,2])

  call pam(x, 2, result)
  write(*,'(a,f10.5)') 'PAM mean dissimilarity: ', result%objective
  write(*,'(a,*(i0,1x))') 'PAM medoids: ', result%medoids
  write(*,'(a,*(i0,1x))') 'PAM clustering: ', result%clustering

  call clara(x, 2, result, samples=5, sample_size=6, seed=123)
  write(*,'(a,f10.5)') 'CLARA mean dissimilarity: ', result%objective

  call fanny(x, 2, result)
  write(*,'(a)') 'FANNY memberships:'
  do i = 1, size(x,1)
    write(*,'(i3,2f10.5)') i, result%membership(i,:)
  end do
end program example_partition
