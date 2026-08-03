! SPDX-License-Identifier: GPL-2.0-or-later
program demo_cluster
  use cluster, only: dp, partition_result, hierarchy_result, silhouette_result, &
    daisy, pam, agnes, silhouette
  implicit none

  real(dp) :: x(10,2)
  real(dp), allocatable :: distances(:, :)
  integer :: status
  character(len=:), allocatable :: message
  type(partition_result) :: partition
  type(hierarchy_result) :: hierarchy
  type(silhouette_result) :: sil

  x = reshape([0.0_dp,0.2_dp,-0.1_dp,0.3_dp,0.1_dp,5.0_dp,5.2_dp,4.9_dp,5.1_dp,4.8_dp, &
               0.0_dp,-0.2_dp,0.1_dp,0.2_dp,-0.1_dp,5.0_dp,4.8_dp,5.2_dp,5.1_dp,4.9_dp], [10,2])

  call pam(x, 2, partition)
  call daisy(x, distances, status=status, message=message)
  call silhouette(partition%clustering, distances, sil)
  call agnes(x, hierarchy, method='average')

  write(*,'(a,*(i0,1x))') 'PAM medoids: ', partition%medoids
  write(*,'(a,f9.5)') 'PAM objective: ', partition%objective
  write(*,'(a,f9.5)') 'Average silhouette: ', sil%average_width
  write(*,'(a,f9.5)') 'AGNES coefficient: ', hierarchy%coefficient
end program demo_cluster
