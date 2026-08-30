! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation derived from R package spdep 1.4-2.
! Translated/modified 2026-08-30 for modern free-form Fortran.
program test_statistics
   use spdep, only : dp, neighbor_list, spatial_weights, spatial_test_result, &
      local_stat_result, eb_result, cell2nb, nb2listw, moran, moran_test, moran_mc, &
      geary, geary_test, geary_mc, local_moran, local_g, losh, lee, joincount_test, &
      ebest, eblocal, choynowski
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none

   type(neighbor_list) :: nb
   type(spatial_weights) :: listw
   type(spatial_test_result) :: test1
   type(spatial_test_result) :: test2
   type(local_stat_result) :: local
   type(eb_result) :: eb
   real(dp) :: x(4)
   real(dp) :: y(4)
   real(dp) :: cases(4)
   real(dp) :: population(4)
   real(dp), allocatable :: probs(:)
   integer :: binary(4)
   integer :: icases(4)
   integer :: failures

   failures = 0
   nb = cell2nb(1, 4)
   listw = nb2listw(nb, "W")
   x = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp]
   y = [4.0_dp, 1.0_dp, 3.0_dp, 2.0_dp]

   call assert_close(moran(x, listw), 0.29130434782608694_dp, 1.0e-12_dp, &
      "Moran I", failures)
   call assert_close(geary(x, listw), 0.3847826086956522_dp, 1.0e-12_dp, &
      "Geary C", failures)

   test1 = moran_test(x, listw, .true.)
   call assert_close(test1%expectation, -1.0_dp / 3.0_dp, 1.0e-12_dp, &
      "Moran expectation", failures)
   call assert_close(test1%variance, 0.1448624238605335_dp, 1.0e-12_dp, &
      "Moran variance", failures)

   test1 = moran_mc(x, listw, 99, 314159)
   test2 = moran_mc(x, listw, 99, 314159)
   call assert_close(test1%p_value, test2%p_value, 0.0_dp, "Moran MC reproducibility", failures)
   call assert_close(test1%expectation, test2%expectation, 0.0_dp, &
      "Moran MC mean reproducibility", failures)

   test1 = geary_test(x, listw, .true.)
   call assert_close(test1%expectation, 1.0_dp, 1.0e-12_dp, "Geary expectation", failures)
   test1 = geary_mc(x, listw, 99, 271828)
   call assert_true(test1%p_value > 0.0_dp .and. test1%p_value <= 1.0_dp, &
      "Geary MC p-value", failures)

   local = local_moran(x, listw)
   call assert_true(all(.not. ieee_is_nan(local%statistic)), "local Moran finite statistics", failures)
   local = local_g(x, listw)
   call assert_true(all(.not. ieee_is_nan(local%statistic)), "local G finite statistics", failures)
   local = losh(x, listw)
   call assert_true(all(.not. ieee_is_nan(local%statistic)), "LOSH finite statistics", failures)
   call assert_true(.not. ieee_is_nan(lee(x, y, listw)), "Lee association finite", failures)

   binary = [1, 1, 0, 0]
   test1 = joincount_test(binary, nb2listw(nb, "B"))
   call assert_close(test1%statistic, 1.0_dp, 1.0e-12_dp, "join count", failures)

   cases = [5.0_dp, 8.0_dp, 12.0_dp, 6.0_dp]
   population = [100.0_dp, 120.0_dp, 180.0_dp, 90.0_dp]
   eb = ebest(cases, population)
   call assert_true(all(eb%estimate >= 0.0_dp), "global EB estimates", failures)
   eb = eblocal(cases, population, nb)
   call assert_true(all(eb%estimate >= 0.0_dp), "local EB estimates", failures)

   icases = [5, 8, 12, 6]
   probs = choynowski(icases, population)
   call assert_true(all(probs >= 0.0_dp .and. probs <= 1.0_dp), &
      "Choynowski probabilities", failures)

   if (failures /= 0) error stop "test_statistics failed"
   print '(a)', "test_statistics: PASS"

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
      real(dp), intent(in) :: expected !! Reference scalar value expected from the deterministic case.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute difference between actual and expected.
      character(len=*), intent(in) :: label !! Human-readable test label printed when the assertion fails.
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented on failure.

      if (abs(actual - expected) > tolerance) then
         failures = failures + 1
         print '(a,1x,a,2(1x,es24.16))', "FAIL:", trim(label), actual, expected
      end if
   end subroutine assert_close

end program test_statistics
