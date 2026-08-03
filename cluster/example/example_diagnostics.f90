! SPDX-License-Identifier: GPL-2.0-or-later
program example_diagnostics
  use cluster, only: dp, silhouette_result, ellipsoid_result, daisy, silhouette, &
    ellipsoidhull, predict_ellipsoid
  implicit none

  real(dp) :: x(8,2)
  real(dp), allocatable :: distances(:, :), d2(:)
  integer :: labels(8), status
  logical, allocatable :: inside(:)
  character(len=:), allocatable :: message
  type(silhouette_result) :: sil
  type(ellipsoid_result) :: ell

  x = reshape([0.0_dp,0.2_dp,-0.1_dp,0.1_dp,5.0_dp,5.2_dp,4.9_dp,5.1_dp, &
               0.0_dp,-0.1_dp,0.2_dp,0.1_dp,5.0_dp,4.8_dp,5.2_dp,5.1_dp], [8,2])
  labels = [1,1,1,1,2,2,2,2]

  call daisy(x, distances, status=status, message=message)
  call silhouette(labels, distances, sil)
  write(*,'(a,f9.5)') 'Average silhouette width: ', sil%average_width

  call ellipsoidhull(x, ell)
  call predict_ellipsoid(ell, x, d2, inside, status)
  write(*,'(a,f12.5)') 'Enclosing ellipsoid volume: ', ell%volume
  write(*,'(a,l1)') 'All observations enclosed: ', all(inside)
end program example_diagnostics
