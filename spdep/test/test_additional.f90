! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation derived from R package spdep 1.4-2.
! Translated/modified 2026-08-30 for modern free-form Fortran.
program test_additional
   use spdep, only : dp, neighbor_list, spatial_weights, cell2nb, nblag, &
      rotation, complement_nb, nblag_cumul, nb2blocknb, nb2listwdist, &
      autocov_dist, local_geary, nb2listw
   implicit none

   type(neighbor_list) :: nb
   type(neighbor_list) :: cmp
   type(neighbor_list) :: cumul
   type(neighbor_list) :: block_nb
   type(neighbor_list) :: expanded
   type(neighbor_list), allocatable :: lags(:)
   type(spatial_weights) :: listw
   real(dp), allocatable :: rotated(:, :)
   real(dp), allocatable :: ac(:)
   real(dp), allocatable :: lc(:)
   real(dp) :: xy(2, 2)
   real(dp) :: coords4(4, 2)
   real(dp) :: coords3(3, 2)
   real(dp) :: z(3)
   integer :: ids(5)
   integer :: failures

   failures = 0

   xy = 0.0_dp
   xy(1, 1) = 1.0_dp
   xy(2, 2) = 1.0_dp
   rotated = rotation(xy, 0.5_dp * acos(-1.0_dp))
   call assert_close(rotated(1, 1), 0.0_dp, 1.0e-12_dp, "rotation x", failures)
   call assert_close(rotated(1, 2), 1.0_dp, 1.0e-12_dp, "rotation y", failures)

   nb = cell2nb(1, 4)
   cmp = complement_nb(nb)
   call assert_true(all(cmp%neighbors(1)%values == [1, 3, 4]), "graph complement", failures)

   lags = nblag(nb, 2)
   cumul = nblag_cumul(lags)
   call assert_true(all(cumul%neighbors(1)%values == [2, 3]), "cumulative graph lags", failures)

   block_nb = cell2nb(1, 3)
   ids = [1, 1, 2, 3, 3]
   expanded = nb2blocknb(block_nb, ids)
   call assert_true(all(expanded%neighbors(1)%values == [2, 3]), "block expansion edge set", failures)
   call assert_true(all(expanded%neighbors(3)%values == [1, 2, 4, 5]), "block expansion center", failures)

   coords4 = 0.0_dp
   coords4(:, 1) = [0.0_dp, 1.0_dp, 3.0_dp, 6.0_dp]
   listw = nb2listwdist(nb, coords4, "idw", "W")
   call assert_close(listw%weights(2)%values(1), 2.0_dp / 3.0_dp, 1.0e-12_dp, &
      "distance weights first", failures)
   call assert_close(listw%weights(2)%values(2), 1.0_dp / 3.0_dp, 1.0e-12_dp, &
      "distance weights second", failures)

   coords3 = 0.0_dp
   coords3(:, 1) = [0.0_dp, 1.0_dp, 3.0_dp]
   z = [1.0_dp, 2.0_dp, 4.0_dp]
   ac = autocov_dist(z, coords3, 2.0_dp, "inverse", "B")
   call assert_close(ac(1), 2.0_dp, 1.0e-12_dp, "autocov first", failures)
   call assert_close(ac(2), 3.0_dp, 1.0e-12_dp, "autocov second", failures)
   call assert_close(ac(3), 1.0_dp, 1.0e-12_dp, "autocov third", failures)

   listw = nb2listw(cell2nb(1, 3), "W")
   lc = local_geary(z, listw)
   call assert_close(lc(1), 3.0_dp / 7.0_dp, 1.0e-12_dp, "Local Geary first", failures)
   call assert_close(lc(2), 15.0_dp / 14.0_dp, 1.0e-12_dp, "Local Geary second", failures)
   call assert_close(lc(3), 12.0_dp / 7.0_dp, 1.0e-12_dp, "Local Geary third", failures)

   if (failures /= 0) error stop "test_additional failed"
   print '(a)', "test_additional: PASS"

contains

   subroutine assert_true(condition, label, failures)
      logical, intent(in) :: condition !! Boolean condition that must be true for the test to pass.
      character(len=*), intent(in) :: label !! Human-readable test label printed when the assertion fails.
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented on failure.

      if (.not. condition) then
         failures = failures + 1
         print '(a,1x,a)', "FAIL:", trim(label)
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, label, failures)
      real(dp), intent(in) :: actual !! Computed scalar value being checked.
      real(dp), intent(in) :: expected !! Deterministic reference value expected by the test.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute difference between actual and expected.
      character(len=*), intent(in) :: label !! Human-readable test label printed when the assertion fails.
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented on failure.

      if (abs(actual - expected) > tolerance) then
         failures = failures + 1
         print '(a,1x,a,2(1x,es24.16))', "FAIL:", trim(label), actual, expected
      end if
   end subroutine assert_close

end program test_additional
