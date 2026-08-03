! SPDX-License-Identifier: GPL-2.0-or-later
program test_hierarchy
  use cluster, only: dp, hierarchy_result, mona_result, agnes, diana, mona
  implicit none

  real(dp) :: x(6,2)
  integer :: xb(6,3)
  type(hierarchy_result) :: h
  type(mona_result) :: m

  x = reshape([0.0_dp, 0.1_dp, -0.2_dp, 5.0_dp, 5.1_dp, 4.9_dp, &
               0.0_dp, -0.2_dp, 0.1_dp, 5.0_dp, 5.2_dp, 4.8_dp], [6,2])
  call agnes(x, h, method='average')
  call check(h%ok(), 'AGNES status')
  call check(size(h%merge,1) == 5 .and. size(h%order) == 6, 'AGNES dimensions')
  call check(h%coefficient >= 0.0_dp .and. h%coefficient <= 1.0_dp, 'AGNES coefficient')
  call check(all(h%height(2:) >= h%height(:4)-1.0e-12_dp), 'AGNES monotone heights')

  call diana(x, h)
  call check(h%ok(), 'DIANA status')
  call check(size(h%merge,1) == 5 .and. size(h%order) == 6, 'DIANA dimensions')
  call check(h%height(5) >= maxval(h%height(:4)), 'DIANA root diameter')
  call check(h%coefficient >= 0.0_dp .and. h%coefficient <= 1.0_dp, 'DIANA coefficient')

  xb = reshape([0,0,0,1,1,1, 0,0,1,0,1,1, 0,1,0,1,0,1], [6,3])
  call mona(xb, m, max_clusters=4)
  call check(m%ok(), 'MONA status')
  call check(m%n_clusters >= 2 .and. m%n_clusters <= 4, 'MONA cluster count')
  call check(size(m%order) == 6, 'MONA order')

  print '(a)', 'test_hierarchy: PASS'

contains
  subroutine check(condition, text)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: text
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//text
      error stop 1
    end if
  end subroutine check
end program test_hierarchy
