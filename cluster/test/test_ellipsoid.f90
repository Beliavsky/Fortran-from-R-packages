! SPDX-License-Identifier: GPL-2.0-or-later
program test_ellipsoid
  use cluster, only: dp, ellipsoid_result, ellipsoidhull, predict_ellipsoid, &
    ellipsoid_points, volume_ellipsoid, cluster_success
  implicit none

  real(dp) :: x(8,2)
  real(dp), allocatable :: distance2(:), boundary(:, :)
  logical, allocatable :: inside(:)
  integer :: status
  type(ellipsoid_result) :: result

  x = reshape([1.0_dp,-1.0_dp,0.0_dp,0.0_dp,0.7_dp,-0.7_dp,0.7_dp,-0.7_dp, &
               0.0_dp,0.0_dp,2.0_dp,-2.0_dp,1.4_dp,-1.4_dp,-1.4_dp,1.4_dp], [8,2])
  call ellipsoidhull(x, result, tolerance=1.0e-6_dp)
  call check(result%ok(), 'ellipsoidhull status')
  call check(volume_ellipsoid(result) > 0.0_dp, 'ellipsoid volume')
  call predict_ellipsoid(result, x, distance2, inside, status)
  call check(status == cluster_success, 'ellipsoid prediction status')
  call check(all(inside), 'input points inside ellipsoid')
  call ellipsoid_points(result, 32, boundary, status)
  call check(status == cluster_success .and. size(boundary,1) == 32, 'ellipsoid points')
  call predict_ellipsoid(result, boundary, distance2, inside, status)
  call check(maxval(abs(distance2-1.0_dp)) < 1.0e-6_dp, 'boundary distances')

  print '(a)', 'test_ellipsoid: PASS'

contains
  subroutine check(condition, text)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: text
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//text
      error stop 1
    end if
  end subroutine check
end program test_ellipsoid
